import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import glob, os, sys
from drugnome_ai.config_class import Config
import textwrap


class BorutaScoreAggregator:

    def __init__(self, cfg, search_dir):
        self.cfg = cfg
        self.search_dir = search_dir
        self.feat_colors = {'Confirmed': '#41AB5D', 'Tentative': '#FFFF00', 'Rejected': '#E31A1C', 'Shadow': '#2171B5'}


    def get_imp_dfs(self):
        self.boruta_imp_df = pd.DataFrame()
        self.total_runs = len(glob.glob(self.search_dir + '/boruta_imp_df*'))
        print('Total runs:', self.total_runs)

        for b_df_file in glob.glob(self.search_dir + '/boruta_imp_df*'):
            #print(b_df_file)

            b_df = pd.read_csv(b_df_file, index_col=False, sep=' ')
            b_df = b_df.reindex(sorted(b_df.columns), axis=1)
            #print(b_df.shape)

            if len(self.boruta_imp_df) > 0:
                self.boruta_imp_df = pd.concat([self.boruta_imp_df, b_df], axis=0)
            else:
                self.boruta_imp_df = b_df
 
        self.boruta_imp_df.to_csv(str(self.cfg.boruta_tables_dir / 'merged.boruta_imp_df.txt'), header=True, index=False)
        print(self.boruta_imp_df.head())
        print(self.boruta_imp_df.columns)
    # --------------------------------------------------------------------------------------------

    def get_final_decisions(self):
   
        if set(['shadowMax', 'shadowMean', 'shadowMin']).issubset(self.boruta_imp_df.columns):
            features = self.boruta_imp_df.drop(['shadowMax', 'shadowMean', 'shadowMin'], axis=1).columns.values
        else:
            features = self.boruta_imp_df.columns.values

        feat_df = pd.DataFrame(0, index=features, columns=['Confirmed', 'Tentative'])

        # read confirmed files
        for confirmed_file in glob.glob(self.search_dir + '/confirmed*'):
            conf_df = pd.read_csv(confirmed_file, sep='\n', header=None, index_col=False)
            conf_series = conf_df[0]

            conf_series = feat_df.index.intersection(conf_series)
            feat_df.loc[conf_series, 'Confirmed'] = feat_df.loc[conf_series, 'Confirmed'] + 1

        # read tentative files
        for tentative_file in glob.glob(self.search_dir + '/tentative*'):
            tent_df = pd.read_csv(tentative_file, sep='\n', header=None, index_col=False)
            tent_series = tent_df[0]

            tent_series = feat_df.index.intersection(tent_series)
            feat_df.loc[tent_series, 'Tentative'] = feat_df.loc[tent_series, 'Tentative'] + 1

        feat_df['Rejected'] = self.total_runs - (feat_df['Confirmed'] + feat_df['Tentative'])
        feat_df = (feat_df / self.total_runs) * 100
        feat_df.sort_values(by=['Confirmed', 'Tentative'], inplace=True)


        fig, ax = plt.subplots(figsize=(20, 10))
        ax = feat_df.plot(kind='bar', stacked=True, figsize=(20, 10), color=[self.feat_colors['Confirmed'], self.feat_colors['Tentative'], self.feat_colors['Rejected']])
        ax.legend(bbox_to_anchor=(1.0, 0.5))

        _ = ax.axhline(self.cfg.boruta_decision_thres, linewidth=0.5, linestyle='--', color='#bdbdbd')
        _ = ax.set_title('Feature Importance Decisions distirbution across all runs')
        _ = ax.set_xlabel('Features')
        _ = ax.set_ylabel('% Decision per Class')
        ax.get_figure().savefig(str(self.cfg.boruta_figs_dir / 'boruta_stacked_barplot.pdf'), bbox_inches='tight')
        # plt.show()


        final_confirmed_features = feat_df.loc[feat_df.Confirmed >= self.cfg.boruta_decision_thres].index.values
        rejected_decision_thres = 90
        final_rejected_features = feat_df.loc[feat_df.Rejected >= rejected_decision_thres].index.values
        final_tentative_features = np.setdiff1d(feat_df.index.values,
                                                np.union1d(final_confirmed_features, final_rejected_features))

        feat_df['Final_Decision'] = 0
        feat_df.loc[final_confirmed_features, 'Final_Decision'] = 'Confirmed'
        feat_df.loc[final_tentative_features, 'Final_Decision'] = 'Tentative'
        feat_df.loc[final_rejected_features, 'Final_Decision'] = 'Rejected'

        self.final_decision_colors = pd.Series([self.feat_colors[feat_df.loc[i, 'Final_Decision']] for i in feat_df.index.values],
                                          index=feat_df.index.values)


        feature_importance = self.boruta_imp_df.mean()

        # Store confirmed/tentative/rejected features into files - to read for feature selection e.g. by Stacking or other classifiers too
        for final_decision in ['Confirmed', 'Tentative', 'Rejected']:
            attributes = feat_df.loc[feat_df.Final_Decision == final_decision, :].index.tolist()
            attributes_importance = feature_importance[attributes].sort_values(ascending=False)
            print("\nattributes_importance")
            print(attributes_importance.head())
            print('\n')
            attributes_importance.to_csv(str(self.cfg.boruta_tables_dir / (final_decision + '.boruta_features.csv')), index=True)



    def get_boruta_boxplots(self):

        final_df = self.boruta_imp_df.reindex(self.boruta_imp_df.median().sort_values().index, axis=1)

        final_decision_colors = self.final_decision_colors.reindex(self.boruta_imp_df.median().sort_values().index)
        final_decision_colors.shadowMin = self.feat_colors['Shadow']
        final_decision_colors.shadowMean = self.feat_colors['Shadow']
        final_decision_colors.shadowMax = self.feat_colors['Shadow']

        for color in final_decision_colors.unique():

            for category, symbol in self.feat_colors.items():
                if symbol == color:
                    feature_cat = category

            final_decision_colors_copy = final_decision_colors[final_decision_colors == color]
            fig, ax = plt.subplots(figsize=(20, 10))
            legend_colors = []
            legend_labels = []

            max_shadowMax = max(final_df['shadowMax'])

            for position, name in enumerate(final_decision_colors_copy.index.values):

                cur_feature_series = final_df[name]
                cur_feature_series = cur_feature_series[~cur_feature_series.isnull()]

                bp = ax.boxplot(cur_feature_series, positions=[position], patch_artist=True, notch=True, widths=0.4,
                                vert=False, flierprops=dict(marker='o', markerfacecolor='black', markersize=3,
                                                linestyle='dotted'))

                cur_face_color = final_decision_colors_copy[name]
                for element in ['boxes', 'whiskers', 'fliers', 'caps']:
                    _ = plt.setp(bp[element], color=cur_face_color)
                _ = plt.setp(bp['medians'], color='#737373')

                for patch in bp['boxes']:
                    _ = patch.set(facecolor=cur_face_color, alpha=0.9)

            labels = [label.replace("_", " ") for label in final_decision_colors_copy.index.values]

            _ = ax.axvline(max_shadowMax, linewidth=0.5, linestyle='--', color=self.feat_colors['Shadow'])
            _ = ax.set_title('Feature importance in disease agnostic run')
            _ = ax.set_yticks(range(position + 1))
            #             _ = ax.set_yticklabels([textwrap.fill(label, 10) for label in labels], fontsize=12, horizontalalignment="center")
            _ = ax.set_yticklabels(labels, fontsize=12)
            _ = ax.tick_params(axis='both', which='major', pad=15)
            #             _ = ax.set_xlim(left=-0.5)
            _ = ax.set_ylabel('Features')
            _ = ax.set_xlabel('Boruta Feature Importance\n(Z-score)')
            plt.tight_layout()
            # plt.show()
            file_name = str(feature_cat) +'_boruta_feature_imp_boxplots.pdf'
            fig.savefig(os.path.join(self.cfg.boruta_figs_dir, file_name), bbox_inches='tight')


if __name__ == '__main__':

    config_file = '../../../config.yaml'
    cfg = Config(config_file)

    agg = BorutaScoreAggregator(cfg, str(cfg.boruta_tables_dir / 'out'))
    # agg = BorutaScoreAggregator(cfg, 'all-mantis-1/feature_selection/boruta/out')
    agg.get_imp_dfs()
    agg.get_final_decisions()
    agg.get_boruta_boxplots()
