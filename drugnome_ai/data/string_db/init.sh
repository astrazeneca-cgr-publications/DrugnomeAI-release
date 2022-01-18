#!/bin/bash
cat drugnome_ai/data/string_db/split_physical_* > drugnome_ai/data/string_db/processed_physical_links.csv.gz
cat drugnome_ai/data/string_db/split_protein_* > drugnome_ai/data/string_db/processed_protein_links.csv.gz
