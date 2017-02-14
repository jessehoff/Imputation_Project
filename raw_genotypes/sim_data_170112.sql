﻿COPY (
SELECT international_id, CASE WHEN sire_international_id IS NULL THEN '0' ELSE sire_international_id END,
	CASE WHEN dam_international_id IS NULL THEN '0' ELSE dam_international_id END,
	CASE WHEN sex = 'M' THEN 1
	WHEN sex = 'F' THEN 2 ELSE 0 END AS sexx,
	array_length(gen_array,1) AS num_loci, gen_array
FROM sample_sheet ss, gen_9913_ggp90kt g
WHERE ss.sample_id = g.sample_id
	AND do_not_analyze IS NULL
	AND international_id LIKE 'SIM%'
	AND assay = 'GGP90KT'
)
TO '/scratch/schnabelr/pg_data_dump/sim_ggp90kt.csv'
WITH CSV DELIMITER '|';
--Query returned successfully: 3498 rows affected, 23.6 secs execution time.
----------------------------------
COPY (
SELECT international_id, CASE WHEN sire_international_id IS NULL THEN '0' ELSE sire_international_id END,
	CASE WHEN dam_international_id IS NULL THEN '0' ELSE dam_international_id END,
	CASE WHEN sex = 'M' THEN 1
	WHEN sex = 'F' THEN 2 ELSE 0 END AS sexx,
	array_length(gen_array,1) AS num_loci, gen_array
FROM sample_sheet ss, gen_9913_ggpf250 g
WHERE ss.sample_id = g.sample_id
	AND do_not_analyze IS NULL
	AND international_id LIKE 'SIM%'
	AND assay = 'GGPF250'
)
TO '/scratch/schnabelr/pg_data_dump/sim_ggpf250.csv'
WITH CSV DELIMITER '|';
--Query returned successfully: 325 rows affected, 6.4 secs execution time.
----------------------------------
COPY (
SELECT international_id, CASE WHEN sire_international_id IS NULL THEN '0' ELSE sire_international_id END,
	CASE WHEN dam_international_id IS NULL THEN '0' ELSE dam_international_id END,
	CASE WHEN sex = 'M' THEN 1
	WHEN sex = 'F' THEN 2 ELSE 0 END AS sexx,
	array_length(gen_array,1) AS num_loci, gen_array
FROM sample_sheet ss, gen_9913_ggphdv3 g
WHERE ss.sample_id = g.sample_id
	AND do_not_analyze IS NULL
	AND international_id LIKE 'SIM%'
	AND assay = 'GGPHDv3'
)
TO '/scratch/schnabelr/pg_data_dump/sim_ggphdv3.csv'
WITH CSV DELIMITER '|';
--Query returned successfully: 2326 rows affected, 25.4 secs execution time.
----------------------------------
COPY (
SELECT international_id, CASE WHEN sire_international_id IS NULL THEN '0' ELSE sire_international_id END,
	CASE WHEN dam_international_id IS NULL THEN '0' ELSE dam_international_id END,
	CASE WHEN sex = 'M' THEN 1
	WHEN sex = 'F' THEN 2 ELSE 0 END AS sexx,
	array_length(gen_array,1) AS num_loci, gen_array
FROM sample_sheet ss, gen_9913_ggpldv3 g
WHERE ss.sample_id = g.sample_id
	AND do_not_analyze IS NULL
	AND international_id LIKE 'SIM%'
	AND assay = 'GGPLDV3'
)
TO '/scratch/schnabelr/pg_data_dump/sim_ggpldv3.csv'
WITH CSV DELIMITER '|';
--Query returned successfully: 3126 rows affected, 6.4 secs execution time.
----------------------------------
COPY (
SELECT international_id, CASE WHEN sire_international_id IS NULL THEN '0' ELSE sire_international_id END,
	CASE WHEN dam_international_id IS NULL THEN '0' ELSE dam_international_id END,
	CASE WHEN sex = 'M' THEN 1
	WHEN sex = 'F' THEN 2 ELSE 0 END AS sexx,
	array_length(gen_array,1) AS num_loci, gen_array
FROM sample_sheet ss, gen_9913_ggpldv4 g
WHERE ss.sample_id = g.sample_id
	AND do_not_analyze IS NULL
	AND international_id LIKE 'SIM%'
	AND assay = 'GGPLDV4'
)
TO '/scratch/schnabelr/pg_data_dump/sim_ggpldv4.csv'
WITH CSV DELIMITER '|';
--Query returned successfully: 2500 rows affected, 6.0 secs execution time.
----------------------------------
COPY (
SELECT international_id, CASE WHEN sire_international_id IS NULL THEN '0' ELSE sire_international_id END,
	CASE WHEN dam_international_id IS NULL THEN '0' ELSE dam_international_id END,
	CASE WHEN sex = 'M' THEN 1
	WHEN sex = 'F' THEN 2 ELSE 0 END AS sexx,
	array_length(gen_array,1) AS num_loci, gen_array
FROM sample_sheet ss, gen_9913_hd g
WHERE ss.sample_id = g.sample_id
	AND do_not_analyze IS NULL
	AND international_id LIKE 'SIM%'
	AND assay = 'HD'
)
TO '/scratch/schnabelr/pg_data_dump/sim_hd.csv'
WITH CSV DELIMITER '|';
--Query returned successfully: 483 rows affected, 34.9 secs execution time.
----------------------------------
COPY (
SELECT international_id, CASE WHEN sire_international_id IS NULL THEN '0' ELSE sire_international_id END,
	CASE WHEN dam_international_id IS NULL THEN '0' ELSE dam_international_id END,
	CASE WHEN sex = 'M' THEN 1
	WHEN sex = 'F' THEN 2 ELSE 0 END AS sexx,
	array_length(gen_array,1) AS num_loci, gen_array
FROM sample_sheet ss, gen_9913_snp50 g
WHERE ss.sample_id = g.sample_id
	AND do_not_analyze IS NULL
	AND international_id LIKE 'SIM%'
	AND assay = 'SNP50'
	AND manifest = 'A'
)
TO '/scratch/schnabelr/pg_data_dump/sim_snp50a.csv'
WITH CSV DELIMITER '|';
--Query returned successfully: 315 rows affected, 1.5 secs execution time.
----------------------------------
COPY (
SELECT international_id, CASE WHEN sire_international_id IS NULL THEN '0' ELSE sire_international_id END,
	CASE WHEN dam_international_id IS NULL THEN '0' ELSE dam_international_id END,
	CASE WHEN sex = 'M' THEN 1
	WHEN sex = 'F' THEN 2 ELSE 0 END AS sexx,
	array_length(gen_array,1) AS num_loci, gen_array
FROM sample_sheet ss, gen_9913_snp50 g
WHERE ss.sample_id = g.sample_id
	AND do_not_analyze IS NULL
	AND international_id LIKE 'SIM%'
	AND assay = 'SNP50'
	AND manifest = 'B'
)
TO '/scratch/schnabelr/pg_data_dump/sim_snp50b.csv'
WITH CSV DELIMITER '|';
--Query returned successfully: 335 rows affected, 1.6 secs execution time.
----------------------------------
COPY (
SELECT international_id, CASE WHEN sire_international_id IS NULL THEN '0' ELSE sire_international_id END,
	CASE WHEN dam_international_id IS NULL THEN '0' ELSE dam_international_id END,
	CASE WHEN sex = 'M' THEN 1
	WHEN sex = 'F' THEN 2 ELSE 0 END AS sexx,
	array_length(gen_array,1) AS num_loci, gen_array
FROM sample_sheet ss, gen_9913_snp50 g
WHERE ss.sample_id = g.sample_id
	AND do_not_analyze IS NULL
	AND international_id LIKE 'SIM%'
	AND assay = 'SNP50'
	AND manifest = 'C'
)
TO '/scratch/schnabelr/pg_data_dump/sim_snp50c.csv'
WITH CSV DELIMITER '|';
--Query returned successfully: 3399 rows affected, 15.6 secs execution time.
----------------------------------

reformat_gen_array_v0.1.1.pl --input sim_ggp90kt.csv --output_prefix GGP90KT 
reformat_gen_array_v0.1.1.pl --input sim_ggpf250.csv --output_prefix GGPF250 
reformat_gen_array_v0.1.1.pl --input sim_ggphdv3.csv --output_prefix GGPHDV3 
reformat_gen_array_v0.1.1.pl --input sim_ggpldv3.csv --output_prefix GGPLDV3 
reformat_gen_array_v0.1.1.pl --input sim_ggpldv4.csv --output_prefix GGPLDV4 
reformat_gen_array_v0.1.1.pl --input sim_snp50a.csv --output_prefix SNP50A 
reformat_gen_array_v0.1.1.pl --input sim_snp50b.csv --output_prefix SNP50B 
reformat_gen_array_v0.1.1.pl --input sim_snp50c.csv --output_prefix SNP50C 









