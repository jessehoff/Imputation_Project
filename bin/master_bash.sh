#!/bin/bash

bash ./bin/allele_call_rate_stats.sh

bash ./bin/allele_call_rate_filter.sh

bash ./bin/individual_call_rate_stats.sh

bash ./bin/individual_call_rate_filter.sh

bash ./bin/hwe_stats.sh

bash ./bin/hwe_filter.sh
