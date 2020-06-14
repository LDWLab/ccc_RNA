#!/usr/bin/env bash

../Score/bin/PhyMeas.py \
	-a ./data/alns/SSU_A-E-noASG.fas \
	-nc \
	-csv \
	-o ./data/cons/SSU_AnoASG-E

./bin/ccc_RNA_command_line.py ./data/alns/SSU_A-E-ASG.fas ./data/cons/SSU_AnoASG-E.csv \
	-o ./output/twc_AnoASG-E/nucs_SSU_A-E-ASG_trunc_by_-2_AnoASG-E \
	-tlc -2 \
	-anno "SSUe_559292_YEAST/1-1800" \
	-traln

./bin/ccc_RNA_command_line.py ./data/alns/SSU_A-E-ASG.fas ./data/cons/SSU_AnoASG-E.csv \
	-o ./output/twc_AnoASG-E/SSU_A-E-ASG_trunc_by_-2_AnoASG-E \
	-tlc -2 \
	-bpfr3d ./data/contacts/SC_LSU_BasePairs_fr3d.csv \
	-chains A2 \
	-anno "SSUe_559292_YEAST/1-1800" \
	-traln

../Score/bin/PhyMeas.py \
	-a ./data/alns/LSU_A-E-noASG.fas \
	-nc \
	-csv \
	-o ./data/cons/LSU_AnoASG-E

./bin/ccc_RNA_command_line.py ./data/alns/LSU_A-E-ASG.fas ./data/cons/LSU_AnoASG-E.csv \
	-o ./output/twc_AnoASG-E/nucs_LSU_A-E-ASG_trunc_by_-2_AnoASG-E \
	-tlc -2 \
	-anno "LSUe_559292_YEAST/1-3554" \
	-traln

./bin/ccc_RNA_command_line.py ./data/alns/LSU_A-E-ASG.fas ./data/cons/LSU_AnoASG-E.csv \
	-o ./output/twc_AnoASG-E/LSU_A-E-ASG_trunc_by_-2_AnoASG-E \
	-tlc -2 \
	-bpfr3d ./data/contacts/SC_LSU_BasePairs_fr3d.csv \
	-chains A4 A1 \
	-anno "LSUe_559292_YEAST/1-3554" \
	-traln

#Uncorrected

../Score/bin/PhyMeas.py \
	-a ./data/alns/LSU_A-E.fas \
	-nc \
	-csv \
	-o ./data/cons/LSU_A+ASG-E

./bin/ccc_RNA_command_line.py ./data/alns/LSU_A-E-ASG.fa ./data/cons/LSU_A+ASG-E.csv \
	-o ./output/uncorrected/nucs_LSU_A-E-ASG_trunc_by_-2_A-E \
	-tlc -2 \
	-anno "LSUe_559292_YEAST/1-3554" \
	-traln

./bin/ccc_RNA_command_line.py ./data/alns/LSU_A-E-ASG.fa ./data/cons/LSU_A+ASG-E.csv \
	-o ./output/uncorrected/LSU_A-E-ASG_trunc_by_-2_A-E \
	-tlc -2 \
	-bpfr3d ./data/contacts/SC_LSU_BasePairs_fr3d.csv \
	-chains A4 A1 \
	-anno "LSUe_559292_YEAST/1-3554" \
	-traln

../Score/bin/PhyMeas.py \
	-a ./data/alns/SSU_A+ASG-E.fas \
	-nc \
	-csv \
	-o ./data/cons/SSU_A+ASG-E

./bin/ccc_RNA_command_line.py ./data/alns/SSU_A-E-ASG.fas ./data/cons/SSU_A+ASG-E.csv \
	-o ./output/uncorrected/nucs_SSU_A-E-ASG_trunc_by_-2_A-E \
	-tlc -2 \
	-anno "SSUe_559292_YEAST/1-1800" \
	-traln

./bin/ccc_RNA_command_line.py ./data/alns/SSU_A-E-ASG.fas ./data/cons/SSU_A+ASG-E.csv \
	-o ./output/uncorrected/SSU_A-E-ASG_trunc_by_-2_A-E \
	-tlc -2 \
	-bpfr3d ./data/contacts/SC_LSU_BasePairs_fr3d.csv \
	-chains A2 \
	-anno "SSUe_559292_YEAST/1-1800" \
	-traln