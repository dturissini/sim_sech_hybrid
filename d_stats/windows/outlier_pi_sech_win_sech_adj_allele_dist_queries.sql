.headers on
.mode tabs

select sum(adj_allele_dist) cum_dist, location, d.sample_id
from outlier_pi_sech_win_sech_adj_allele_dist_50000 d, sample_pop s
where d.sample_id = s.sample_id
group by d.sample_id, location
order by sum(adj_allele_dist);
cum_dist|location|sample_id
0.0|Anro, Seychelles|SECH_3-sech_Anro_B3_TTAGGC_L001
14.0418488294593|Anro, Seychelles|SECH_6-sech_Anro_B8_GCCAAT_L001
15.5623648824455|Anro, Seychelles|SECH_7-sech_Anro_B7_CAGATC_L001
15.759722068179|Anro, Seychelles|SECH_5-sech_Anro_B6_ACAGTG_L001
16.831559487235|La Digue, Seychelles|SECH_14-sech_LD8_AGTTCC_L002
17.1959841898944|Denis, Seychelles|SECH_DENIS_72_ATCGTTAG-TGCGAAGT_S141_L003
20.0283696999954|Marianne, Seychelles|SECH_1-sech_maria_3_ATCACG_L001
21.006943091244|La Digue, Seychelles|SECH_16-sech_LD12_CCGTCC_L002
26.4685798937997|La Digue, Seychelles|SECH_13-sech_LD16_AGTCAA_L002
26.9177658047503|Anro, Seychelles|SECH_9-sech_Anro_B2_GATCAG_L001
28.3963552829538|La Digue, Seychelles|SECH_15-sech_LD13_ATGTCA_L002
31.9175775352629|La Digue, Seychelles|SECH_11-sech_LD15_GGCTAC_L002
32.1100264831236|Anro, Seychelles|SECH_4-sech_Anro_B5_TGACCA_L001
33.2380616685701|La Digue, Seychelles|SECH_10-sech_LD14_TAGCTT_L001
33.5598184138534|Anro, Seychelles|SECH_8-sech_Anro_B1_ACTTGA_L001
38.2656133010898|Marianne, Seychelles|SECH_2-sech_mariane_1_CGATGT_L001
39.395261481793|Anro, Seychelles|SECH_ANRO_71_GTTATTAT-CGCTAAGC_S140_L003
47.9788055621308|La Digue, Seychelles|SECH_12-sech_LD11_CTTGTA_L002
60.720193843125|Denis, Seychelles|DenisNoni10
69.9341478808783|Denis, Seychelles|DenisNF123_1_sampleo3
104.465209081045|Denis, Seychelles|DenisNF100
2661.59108836679|Unknown|14021-0428.25_SRR869587
2980.66769061977|Denis, Seychelles|DenisAT3
3367.39846703708|Denis, Seychelles|DenisDNJ6
3385.83842193727|Denis, Seychelles|DenisNF155_sampleo14
3475.87616082785|Denis, Seychelles|DenisNoni101
3575.76741729128|Denis, Seychelles|DenisNF134
3578.65718133003|Denis, Seychelles|Denis135_sampleo5
3580.2744672527|Denis, Seychelles|DenisMCL_sampleo15
3581.40572811332|Denis, Seychelles|DenisJT1_sampleo13
3593.41645822523|Denis, Seychelles|DenisNoni60
3610.61118598348|Denis, Seychelles|DenisNF66_sampleo20
3621.57936969556|Denis, Seychelles|DenisNF13
3662.67258244888|Praslin, Seychelles|SECH_19-sech_PNF4_GTGAAA_L002
3667.76276902591|Denis, Seychelles|DenisNF13_2
3685.58173700433|Praslin, Seychelles|SECH_23-sech_PNF11_GAGTGG_L001
3715.06860266193|Denis, Seychelles|Denis7_2
3715.06860266193|Denis, Seychelles|Denis7_8_sampleo11
3755.59630585266|Praslin, Seychelles|SECH_21-sech_PNF3_GTTTCG_L001
3758.84586244232|Praslin, Seychelles|SECH_17-sech_PNF8_ATTCCT_L002
3785.13318948686|Denis, Seychelles|Denis124
3822.78651662517|Denis, Seychelles|DenisAMT_sampleo22
3881.06601036041|Praslin, Seychelles|SECH_20-sech_PNF5_GTGGCC_L001
3888.57827240727|Praslin, Seychelles|SECH_22-sech_PNF10_CGTACG_L001
3899.03874574876|Denis, Seychelles|DenisMCL_sampleo18
4078.74726194899|Praslin, Seychelles|SECH_18-sech_PNF7_GTCCGC_L002


select sum(max_dist)
from (select dsw_id, max(adj_allele_dist) max_dist
      from outlier_pi_sech_win_sech_adj_allele_dist_50000 d
      group by dsw_id) x;
sum(max_dist)
6320.60323307753


