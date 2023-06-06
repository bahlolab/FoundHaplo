#!/usr/bin/env bash

FoundHaplo_DIR=$1
TEST_SAMPLES_FILE=$2  # .txt file with test sample IDs # example: FoundHaplo/input_files/input_vcf_data/test_cohort/samples.txt
CONTROL_SAMPLES_FILE=$3 # .txt file with control sample IDs # example: FoundHaplo/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_samples_by_population/EUR.txt

CHUNK_SIZE=$4 # Number of samples in one chunk recomended a max of 1000.

echo "splitting samples files to chunks of "$CHUNK_SIZE" samples"

TEST_SAMPLES_DIR="$(dirname "$TEST_SAMPLES_FILE")"
CONTROL_SAMPLES_DIR="$(dirname "$CONTROL_SAMPLES_FILE")"

mkdir -p $TEST_SAMPLES_DIR/samples
mkdir -p $CONTROL_SAMPLES_DIR/samples

split -l $CHUNK_SIZE -d --additional-suffix=.txt $TEST_SAMPLES_FILE  $TEST_SAMPLES_DIR/samples/file

echo "Control samples will be in chunks of 100 sample ids by default "
split -l 100 -d --additional-suffix=.txt $CONTROL_SAMPLES_FILE  $CONTROL_SAMPLES_DIR/samples/file

echo "Sample IDs for the test cohort are in $TEST_SAMPLES_DIR/samples"
echo "Sample IDs for the control cohort are in $CONTROL_SAMPLES_DIR/samples"


