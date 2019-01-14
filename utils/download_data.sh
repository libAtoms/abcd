#!/bin/bash
set -euo pipefail
mkdir -p data; cd data

for i in \
  TungstenTrainingConfigurations/GAP_{1..6}.xyz \
  IronGAP/Fe_bcc_bulk_vac_multivac_surf_gamma_int_diint.xyz.gz \
  IronGAP/bcc_primitive_high.xyz \
  IronGAP/bcc_primitive_expanded_high.xyz \
  IronGAP/bcc_primitive_contracted_high.xyz \
  IronGAP/bcc_bulk_54_high.xyz \
  IronGAP/bcc_bulk_54_expanded_high.xyz \
  IronGAP/bcc_bulk_54_expanded_2_high.xyz \
  IronGAP/bcc_bulk_128_high.xyz \
  IronGAP/bcc_bulk_128_expanded_high.xyz \
  IronGAP/bcc_monovacancy_53_high.xyz \
  IronGAP/bcc_monovacancy_127_high.xyz \
  IronGAP/bcc_doublevacancy_126_high.xyz \
  IronGAP/bcc_doublevacancy_126_2NN_high.xyz \
  IronGAP/bcc_doublevacancy_126_1NN_high.xyz \
  IronGAP/bcc_trivacancy_100_125_high.xyz \
  IronGAP/bcc_trivacancy_110_125_high.xyz \
  IronGAP/bcc_trivacancy_111_125_high.xyz \
  IronGAP/bcc_quadvacancy_124_high.xyz \
  IronGAP/bcc_quinvacancy_123_high.xyz \
  IronGAP/bcc_self_interstitial_high.xyz \
  IronGAP/bcc_self_di_interstitial_npc_130_high.xyz \
  IronGAP/bcc_surface.xyz \
  IronGAP/bcc_gamma_surface.xyz
do
  wget -nv -nc http://www.libatoms.org/pub/Home/$i
done

for i in *.gz
do
  gunzip -f -k $i
done
