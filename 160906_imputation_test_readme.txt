COPY (
WITH a AS (
	SELECT lab_id, COUNT(assay) AS num_assays
	FROM sample_sheet
	WHERE assay IN ('SNP50','HD','GGPF250')
		AND breed = 100
		AND international_id LIKE 'AAN%'
	GROUP BY lab_id
)
SELECT international_id, CASE WHEN sire_international_id IS NULL THEN '0' ELSE sire_international_id END,
	CASE WHEN dam_international_id IS NULL THEN '0' ELSE dam_international_id END, sex, array_length(gen_array,1), gen_array
FROM a, sample_sheet ss, gen_9913_snp50 g 
WHERE a.lab_id = ss.lab_id
	AND ss.sample_id = g.sample_id
	AND ss.assay = 'SNP50'
	AND num_assays = 3
)
TO '/scratch/schnabelr/pg_data_dump/imp_test_snp50.csv'
WITH CSV DELIMITER '|';
--Query returned successfully: 75 rows affected, 1.0 secs execution time.
----------------------------
COPY (
WITH a AS (
	SELECT lab_id, COUNT(assay) AS num_assays
	FROM sample_sheet
	WHERE assay IN ('SNP50','HD','GGPF250')
		AND breed = 100
		AND international_id LIKE 'AAN%'
	GROUP BY lab_id
)
SELECT international_id, CASE WHEN sire_international_id IS NULL THEN '0' ELSE sire_international_id END,
	CASE WHEN dam_international_id IS NULL THEN '0' ELSE dam_international_id END, sex, array_length(gen_array,1), gen_array
FROM a, sample_sheet ss, gen_9913_hd g 
WHERE a.lab_id = ss.lab_id
	AND ss.sample_id = g.sample_id
	AND ss.assay = 'HD'
	AND num_assays = 3
)
TO '/scratch/schnabelr/pg_data_dump/imp_test_hd.csv'
WITH CSV DELIMITER '|';
--Query returned successfully: 75 rows affected, 11.8 secs execution time.

----------------------------
COPY (
WITH a AS (
	SELECT lab_id, COUNT(assay) AS num_assays
	FROM sample_sheet
	WHERE assay IN ('SNP50','HD','GGPF250')
		AND breed = 100
		AND international_id LIKE 'AAN%'
	GROUP BY lab_id
)
SELECT international_id, CASE WHEN sire_international_id IS NULL THEN '0' ELSE sire_international_id END,
	CASE WHEN dam_international_id IS NULL THEN '0' ELSE dam_international_id END, sex, array_length(gen_array,1), gen_array
FROM a, sample_sheet ss, gen_9913_ggpf250 g 
WHERE a.lab_id = ss.lab_id
	AND ss.sample_id = g.sample_id
	AND ss.assay = 'GGPF250'
	AND num_assays = 3
)
TO '/scratch/schnabelr/pg_data_dump/imp_test_ggpf250.csv'
WITH CSV DELIMITER '|';
--Query returned successfully: 75 rows affected, 4.2 secs execution time.

reformat_gen_array_v0.1.1.pl --input imp_test_snp50.csv --output_prefix imp_test 
reformat_gen_array_v0.1.1.pl --input imp_test_hd.csv --output_prefix imp_test 
reformat_gen_array_v0.1.1.pl --input imp_test_ggpf250.csv --output_prefix imp_test 

COPY (
SELECT international_id, CASE WHEN sire_international_id IS NULL THEN '0' ELSE sire_international_id END,
	CASE WHEN dam_international_id IS NULL THEN '0' ELSE dam_international_id END, sex, array_length(gen_array,1), gen_array
FROM sample_sheet ss, gen_9913_snp50 g 
WHERE ss.sample_id = g.sample_id
	AND assay = 'SNP50'
	AND manifest = 'A'
	AND breed = 100
	LIMIT 100
)
TO '/scratch/schnabelr/pg_data_dump/imp_test_snp50_A.csv'
WITH CSV DELIMITER '|';
--Query returned successfully: 100 rows affected, 1.1 secs execution time.

COPY (
SELECT international_id, CASE WHEN sire_international_id IS NULL THEN '0' ELSE sire_international_id END,
	CASE WHEN dam_international_id IS NULL THEN '0' ELSE dam_international_id END, sex, array_length(gen_array,1), gen_array
FROM sample_sheet ss, gen_9913_snp50 g 
WHERE ss.sample_id = g.sample_id
	AND assay = 'SNP50'
	AND manifest = 'B'
	AND breed = 100
	LIMIT 100
)
TO '/scratch/schnabelr/pg_data_dump/imp_test_snp50_B.csv'
WITH CSV DELIMITER '|';
--Query returned successfully: 100 rows affected, 907 msec execution time.

COPY (
SELECT international_id, CASE WHEN sire_international_id IS NULL THEN '0' ELSE sire_international_id END,
	CASE WHEN dam_international_id IS NULL THEN '0' ELSE dam_international_id END, sex, array_length(gen_array,1), gen_array
FROM sample_sheet ss, gen_9913_snp50 g 
WHERE ss.sample_id = g.sample_id
	AND assay = 'SNP50'
	AND manifest = 'C'
	AND breed = 100
	LIMIT 100
)
TO '/scratch/schnabelr/pg_data_dump/imp_test_snp50_C.csv'
WITH CSV DELIMITER '|';
--Query returned successfully: 100 rows affected, 1.2 secs execution time.

COPY (
SELECT international_id, CASE WHEN sire_international_id IS NULL THEN '0' ELSE sire_international_id END,
	CASE WHEN dam_international_id IS NULL THEN '0' ELSE dam_international_id END, sex, array_length(gen_array,1), gen_array
FROM sample_sheet ss, gen_9913_hd g 
WHERE ss.sample_id = g.sample_id
	AND assay = 'HD'
	AND manifest = 'A'
	AND breed = 100
	LIMIT 100
)
TO '/scratch/schnabelr/pg_data_dump/imp_test_hd_A.csv'
WITH CSV DELIMITER '|';
--Query returned successfully: 100 rows affected, 11.6 secs execution time.

COPY (
SELECT international_id, CASE WHEN sire_international_id IS NULL THEN '0' ELSE sire_international_id END,
	CASE WHEN dam_international_id IS NULL THEN '0' ELSE dam_international_id END, sex, array_length(gen_array,1), gen_array
FROM sample_sheet ss, gen_9913_hd g 
WHERE ss.sample_id = g.sample_id
	AND assay = 'HD'
	AND manifest = 'B'
	AND breed = 100
	LIMIT 100
)
TO '/scratch/schnabelr/pg_data_dump/imp_test_hd_B.csv'
WITH CSV DELIMITER '|';
--Query returned successfully: 100 rows affected, 14.1 secs execution time.

COPY (
SELECT international_id, CASE WHEN sire_international_id IS NULL THEN '0' ELSE sire_international_id END,
	CASE WHEN dam_international_id IS NULL THEN '0' ELSE dam_international_id END, sex, array_length(gen_array,1), gen_array
FROM sample_sheet ss, gen_9913_ggpf250 g 
WHERE ss.sample_id = g.sample_id
	AND assay = 'GGPF250'
	AND manifest = 'A'
	AND breed = 100
	LIMIT 100
)
TO '/scratch/schnabelr/pg_data_dump/imp_test_ggpf250_A.csv'
WITH CSV DELIMITER '|';
--Query returned successfully: 100 rows affected, 3.5 secs execution time.

reformat_gen_array_v0.1.1.pl --input imp_test_snp50_A.csv --output_prefix test_snp50_A
reformat_gen_array_v0.1.1.pl --input imp_test_snp50_B.csv --output_prefix test_snp50_B
reformat_gen_array_v0.1.1.pl --input imp_test_snp50_C.csv --output_prefix test_snp50_C
reformat_gen_array_v0.1.1.pl --input imp_test_hd_A.csv --output_prefix test_hd_A
reformat_gen_array_v0.1.1.pl --input imp_test_ggpf250_A.csv --output_prefix test_ggpf250_A


