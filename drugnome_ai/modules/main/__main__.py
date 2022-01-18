import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from argparse import ArgumentParser
import sys, os
import glob
import pandas as pd
import ntpath
import pickle
from argparse import RawTextHelpFormatter
from pathlib import Path

class DrugnomeAI:

	def __init__(self, output_dir, config_file=None, nthreads=4, iterations=10, custom_known_genes_file=None, fast_run_option=False, superv_models=None,
                     data_source=None, tier_tag=None, pharos_tag=None, inter_pro=None, mantis_data_option=None, genic_intol_option=None):

		from drugnome_ai.config_class import Config
		self.config_file = config_file
		self.output_dir = output_dir

		self.cfg = Config(self.output_dir, self.config_file)

		# modify default config paramters when provided with respective parameters
		self.cfg.nthreads = int(nthreads)
		self.cfg.iterations = int(iterations)

		if fast_run_option:
			self.cfg.classifiers = ['ExtraTreesClassifier', 'RandomForestClassifier', 'SVC', 'GradientBoostingClassifier']

		if superv_models:
			models_dict = { 'et': 'ExtraTreesClassifier',
							'rf': 'RandomForestClassifier',
							'svc': 'SVC',
							'gb': 'GradientBoostingClassifier',
							'xgb': 'XGBoost',
							'dnn': 'DNN',
							'stack': 'Stacking',
							'nb': 'GaussianNB'
							}

			try:
				self.cfg.classifiers = list(set([ models_dict[k] for k in superv_models ]))
			except:
				print('[Warning] -m option args are not correct.\n\t  Currently going ahead with drugnome_ai run using the 8 default classifiers (unless -f has also been specified which will integrate 4 classifiers only).\n')


		if data_source:
			data_dict = {
						'pharos': 'pharos',
						'inter': 'inter_pro'
						}
			try:
				self.cfg.data_source = list(set([ data_dict[k] for k in data_source ]))
				print(self.cfg.data_source)
			except:
				self.cfg.data_source = ['pharos', 'inter_pro']
				print('All data sources selected: ', self.cfg.data_source, '\n')

		if 'inter_pro' in self.cfg.data_source:
			pro_dict = {'dom': 'domains',
						'fam': 'families',
						'sup': 'super_families'
						}
			try:
				self.cfg.inter_pro = list(set([ pro_dict[k] for k in inter_pro ]))
				print(self.cfg.inter_pro)
			except:
				self.cfg.inter_pro = ['domains', 'families', 'super_families']
				print('[Warning] -int option args if specified may not be correct.\n\t  Currently going ahead with drugnome-ai running all InterPro domains.\n')


		if custom_known_genes_file:
			self.cfg.custom_known_genes_file = Path(custom_known_genes_file)

		elif tier_tag:
			tier_dict =  { '1': 'Tier 1',
						   '2': 'Tier 2',
						   '3A': 'Tier 3A',
						   '3B': 'Tier 3B'
						   }
			try:
				self.cfg.tier_tag = list(set([ tier_dict[k] for k in tier_tag ]))
				print(self.cfg.tier_tag)
			except:
				self.cfg.tier_tag = ['Tier 1', 'Tier 2', 'Tier 3A', 'Tier 3B']
				print('[Warning] -t option args are not correct.\n\t  Currently going ahead with drugnome-ai run using all 3 tiers.\n')

		elif pharos_tag:
			pharos_dict = { 'tclin': 'Tclin',
                            'tchem' : 'Tchem',
                            'tbio' : 'Tbio',
                            'tdark' : 'Tdark'
                        	}
			try:
				self.cfg.pharos_tag = list(set([ pharos_dict[k] for k in pharos_tag ]))
			except:
				self.cfg.pharos_tag = ['Tclin', 'Tchem', 'Tbio', 'Tdark']
				print('[Warning] -p option args are not correct. Running default all')

		if mantis_data_option:
				self.cfg.mantis_data_option = mantis_data_option

		print('nthreads:', self.cfg.nthreads)
		print('Stochastic iterations:', self.cfg.iterations)
		print('Classifiers:', self.cfg.classifiers)
		print('Custom known genes:', self.cfg.custom_known_genes_file)
		print("Data sources: ", self.cfg.data_source)
		if 'inter_pro' in self.cfg.data_source:
			print("InterPro features: ", self.cfg.inter_pro)

		# Run profiler and store results to ouput dir
		# if self.cfg.disease_oriented:
		# 	os.system("mantisml-profiler -vc " + config_file + " -o " + self.output_dir + " > " + str(self.cfg.out_root) + "/profiler_metadata.out")


	def get_clf_id_with_top_auc(self):

		auc_per_clf = {}

		metric_files = glob.glob(str(self.cfg.superv_out / 'PU_*.evaluation_metrics.tsv'))

		for f in metric_files:
			clf_id = ntpath.basename(f).split('.')[0].replace('PU_', '')

			tmp_df = pd.read_csv(f, sep='\t', index_col=0)
			avg_auc = tmp_df.AUC.median()
			auc_per_clf[clf_id] = avg_auc

		top_clf = max(auc_per_clf, key=auc_per_clf.get)
		print('Top classifier:', top_clf)

		return top_clf

	def run(self, clf_id=None, final_level_classifier='DNN', run_eda=False, run_pu=False,
			run_aggregate_results=False, run_merge_results=False, run_boruta=False,
			run_unsupervised=False):

		# *** Load required modules ***
		from drugnome_ai.modules.supervised_learn.pu_learn.pu_learning import PULearning
		from drugnome_ai.modules.pre_processing.eda_wrapper import EDAWrapper
		from drugnome_ai.modules.pre_processing.feature_table_compiler import FeatureTableCompiler
		from drugnome_ai.modules.pre_processing.agnostic_table_compiler import AgnosticTableCompiler
		from drugnome_ai.modules.unsupervised_learn.dimens_reduction_wrapper import DimensReductionWrapper
		from drugnome_ai.modules.post_processing.process_classifier_results import ProcessClassifierResults
		from drugnome_ai.modules.post_processing.merge_predictions_from_classifiers import MergePredictionsFromClassifiers
		from drugnome_ai.modules.supervised_learn.feature_selection.run_boruta import BorutaWrapper

		# ========= Run EDA and pre-processing =========
		if run_eda:
			# ========= Compile feature table =========
			if not self.cfg.disease_oriented:
				agnostic_compiler = AgnosticTableCompiler(self.cfg)
				agnostic_compiler.run()

			else:
				feat_compiler = FeatureTableCompiler(self.cfg)
				feat_compiler.run()

			eda_wrapper = EDAWrapper(self.cfg)
			eda_wrapper.run()

		data = pd.read_csv(self.cfg.processed_data_dir / "processed_feature_table.tsv", sep='\t')

		# ================== Supervised methods ==================
		# ************ Run PU Learning ************
		if run_pu:
			for clf_id in self.cfg.classifiers:
				print('Classifier:', clf_id)
				pu = PULearning(self.cfg, data, clf_id, final_level_classifier)
				pu.run()

		# ************ Process predictions per classifier ************
		if run_aggregate_results:
			aggr_res = ProcessClassifierResults(self.cfg, show_plots=True)
			aggr_res.run()

		# ************ Merge results from all classifiers ************
		if run_merge_results:
			merger = MergePredictionsFromClassifiers(self.cfg)
			merger.run()

		# ************ Run Boruta feature seleciton algorithm ************
		if run_boruta:
			boru_wrapper = BorutaWrapper(self.cfg)
			boru_wrapper.run()

		# ========= Unsupervised methods =========
		# PCA, sparse PCA and t-SNE
		if run_unsupervised:
			recalc = False # default: False

			if clf_id is None:
				highlighted_genes = self.cfg.highlighted_genes
			else:
				top_genes_num = 40
				novel_genes = pd.read_csv(str(self.cfg.superv_ranked_by_proba / (clf_id + '.Novel_genes.Ranked_by_prediction_proba.csv')), header=None, index_col=0)
				highlighted_genes = novel_genes.head(top_genes_num).index.values

			dim_reduct_wrapper = DimensReductionWrapper(self.cfg, data, highlighted_genes, recalc)
			dim_reduct_wrapper.run()

	def run_non_clf_specific_analysis(self):
		""" run_tag: pre """
		args_dict = {'run_eda': True, 'run_unsupervised': self.cfg.run_unsupervised}
		self.run(**args_dict)

	def run_boruta_algorithm(self):
		""" run_tag: boruta """
		args_dict = {'run_boruta': True}
		self.run(**args_dict)

	def run_pu_learning(self):
		""" run_tag: pu """
		args_dict = {'run_pu': True}
		self.run(**args_dict)

	def run_post_processing_analysis(self):
		""" run_tag: post """
		args_dict = {'run_aggregate_results': True, 'run_merge_results': True}
		self.run(**args_dict)

	def run_clf_specific_unsupervised_analysis(self, clf_id):
		""" run_tag: post_unsup """
		args_dict = {'clf_id': clf_id, 'run_unsupervised': True}
		self.run(**args_dict)

	def run_debug_mode_analysis(self):
		""" run_tag: debug """
		args_dict = {'run_eda': True, 'run_pu': True, 'run_aggregate_results': False, 'run_merge_results': False,
				  'run_boruta': False, 'run_unsupervised': False}
		self.cfg.debug_mode = True
		self.cfg.nthreads = 2
		self.cfg.iterations = 1
		self.cfg.classifiers = ['ExtraTreesClassifier', 'GradientBoostingClassifier']
		self.run(**args_dict)

	# ---------------------- Run Full pipeline ------------------------
	def run_all(self):
		""" run_tag: all """
		args_dict = {'run_eda': True, 'run_pu': True, 'run_aggregate_results': True, 'run_merge_results': True,
				  'run_boruta': False, 'run_unsupervised': True}
		self.run(**args_dict)
	# -----------------------------------------------------------------

def main():
	parser = ArgumentParser(formatter_class=RawTextHelpFormatter)
	parser.add_argument("-o", "--output-dir", dest="output_dir", help="Output directory name\n(absolute/relative path e.g. ./CKD, /tmp/Epilepsy-testing, etc.)\nIf it doesn't exist it will automatically be created [Required]\n\n", required=True)
	parser.add_argument("-c", "--config-file", dest="config_file", help="Config file (.yaml) with run parameters [Required for disease specific analysis]\n\n", default=None)
	parser.add_argument("-r", "--run-tag", dest="run_tag", choices=['all', 'pre', 'boruta', 'pu', 'post', 'post_unsup', 'debug'], default='all', help="Specify type of analysis to run (default: all)\n\n")
	parser.add_argument("-f", "--fast", action="count", help="Fast training using only 4 classifiers: Extra Trees, Random Forest, SVC and Gradient Boosting.\nBy default, mantis-ml uses 6 supervised models for training: Extra Trees, Random Forest, SVC, Gradient Boosting, XGBoost and Deep Neural Net.\n\n")
	parser.add_argument("-s", "--superv-models", dest="superv_models", nargs="+", type=str, choices=['et', 'rf', 'svc', 'gb', 'xgb', 'dnn', 'stack', 'nb'], default=None, help="Explicitly specify which supervised models to be used for training. This overrides the '-f/--fast' option.\n- Options:\n et: Extra Trees\n rf: Random Forest\n gb: Gradient Boosting\n xgb: XGBoost\n svc: Support Vector Classifier\n dnn: Deep Neural Net\n stack: Stacking classifier\n nb: Naive Bayes\n\nMultiple models may be specified using a ',' separator, e.g. -m et,rf,stack\nWhen this option is not specified, 6 models are trained by default with mantis-ml: Extra Trees, Random Forest, SVC, Gradient Boosting, XGBoost and Deep Neural Net. \n\n")
	parser.add_argument("-d", "--data-source", dest="data_source", choices=['pharos', 'inter'], nargs="+", type=str, default=None, help="Specify which data sources to include - pharos: pharos, inter: inter_pro] (default: all)\n\n")
	parser.add_argument("-x", "--inter-pro", dest="inter_pro", nargs="+", type=str, choices=['dom', 'fam', 'sup'], default=None, help="Specify which InterPro features you would like - dom: domains, fam: families, sup: super homologous families. (default: None)\n Labels can be specified either by selecting -t/-p option or by a list of custom seed genes using -k option\n\n")
	parser.add_argument("-t", "--tier-tag", dest="tier_tag", nargs="+", type=str, choices=['1', '2', '3A', '3B'], default=None, help="Specify the tier(s) of druggability you wish to train on. 1: Tier 1, 2: Tier 2, 3A: Tier 3A, and 3B Tier 3B (default: all)\n\n")
	parser.add_argument("-p", "--pharos-tag", dest="pharos_tag", nargs="+", type=str, choices=['tclin', 'tchem', 'tbio', 'tdark'], default=None, help="Specify the pharos tier(s) of druggability you wish to train on. tclin: Tclin, tchem: Tchem, tbio: Tbio, and tdark: Tdark (default: all)\n\n")
	parser.add_argument("-m", "--mantis-ml", dest='mantis_features', action="count", help="Choose to use generic feautres derived from mantis-ml; incl. ExAC, Essential Mouse Genes, GnomAD, Genic Intolerance Scores, GWAS & MGI Essential features.\n\n")
	parser.add_argument("-l", "--genic-intol", dest='genic_intol', action="count", help="Choose to include Genic Intolerance Scores from mantis-ml\n\n")
	parser.add_argument("-k", "--known-genes-file", dest="known_genes_file", help="File with custom list of known genes used for training (new-line separated)\n\n")
	parser.add_argument("-n", "--nthreads", dest="nthreads", default=4, help="Number of threads (default: 4)\n\n")
	parser.add_argument("-i", "--iterations", dest="iterations", default=10, help="Number of stochastic iterations for semi-supervised learning (default: 10)\n\n")

	if len(sys.argv)==1:
		parser.print_help(sys.stderr)
		sys.exit(1)

	args = parser.parse_args()
	print(args)

	if args.tier_tag and args.pharos_tag:
		print('[Warning] - both -p and -t cannot be called together. Program exiting...')
		sys.exit(1)

	if (args.tier_tag or args.pharos_tag) and args.known_genes_file:
		print('[Warning] - both -p/-t and -k cannot be called together. Select a label option or provide your own known genes list. Program exiting...')
		sys.exit(1)

	output_dir = args.output_dir
	config_file = args.config_file
	run_tag = args.run_tag
	fast_run_option = bool(args.fast)
	superv_models = args.superv_models
	data_source = args.data_source
	tier_tag = args.tier_tag
	pharos_tag = args.pharos_tag
	inter_pro = args.inter_pro
	mantis_data_option = bool(args.mantis_features)
	genic_intol_option = bool(args.genic_intol)

	custom_known_genes_file = args.known_genes_file
	nthreads = args.nthreads
	iterations = args.iterations

	drugnome = DrugnomeAI(output_dir,
			  config_file,
			  nthreads=nthreads,
			  iterations=iterations,
			  custom_known_genes_file=custom_known_genes_file,
			  fast_run_option = fast_run_option,
			  superv_models = superv_models,
			  data_source = data_source,
              tier_tag=tier_tag,
              pharos_tag=pharos_tag,
			  mantis_data_option=mantis_data_option,
			  genic_intol_option=genic_intol_option,
              inter_pro = inter_pro)

	if run_tag == 'all':
		drugnome.run_all()
	elif run_tag == 'pre':
		drugnome.run_non_clf_specific_analysis()
	elif run_tag == 'pu':
		drugnome.run_pu_learning()
	elif run_tag == 'post':
		drugnome.run_post_processing_analysis()
	elif run_tag == 'post_unsup':
		top_clf = drugnome.get_clf_id_with_top_auc()
		drugnome.run_clf_specific_unsupervised_analysis(top_clf)
	elif run_tag == 'boruta':
		drugnome.run_boruta_algorithm()
	elif run_tag == 'debug':
		drugnome.run_debug_mode_analysis()

if __name__ == '__main__':
	main()


