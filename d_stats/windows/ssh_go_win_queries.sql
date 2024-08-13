sqlite3 /work/users/d/t/dturissi/drosophila/ssh/introgression/d_stats/windows/ssh_d_win.db

attach database '/proj/matutelb/projects/gwas/genotype_datasets/sech_oa/sech_oa_anno.db' as a;
attach database '/proj/matutelb/projects/gwas/orthodb/odb11v0_subsetted_sim_mel.db' as o;
.headers on
.separator "\t"






select count(distinct gi.ncbi_gene_id) num_genes, 1.0 * count(distinct gi.ncbi_gene_id) / total_genes_windows prop_genes_win, prop_genes_genome, 
       1.0 * count(distinct gi.ncbi_gene_id) / total_genes_windows / prop_genes_genome enrichment, count(distinct pw_id) num_wins, name
from poly_win_50000 w, gene_info gi, sim_go sg, go_edges_tier_2 t, go,
     (select t2.go_id_parent go_id, 1.0 * count(distinct ncbi_gene_id) / total_genes prop_genes_genome
      from sim_go sg2, go_edges_tier_2 t2, (select count(distinct ncbi_gene_id) total_genes from sim_go) x2
      where sg2.go_id = t2.go_id_child
      group by t2.go_id_parent, total_genes) x,
      (select count(distinct ncbi_gene_id) total_genes_windows
       from poly_win_50000 w2, gene_info gi2
       where gi2.start between w2.start and w2.end
       and w2.chrom = gi2.chrom
       and pi_sech > 0.025) p
where pi_sech > 0.025
and gi.start between w.start and w.end
and w.chrom = gi.chrom
and gi.ncbi_gene_id = sg.ncbi_gene_id
and sg.go_id = t.go_id_child
and t.go_id_parent = go.go_id
and namespace = 'biological_process'
and x.go_id = go.go_id
group by name, total_genes_windows, prop_genes_genome
order by enrichment desc, num_genes desc, num_wins desc, name;

num_genes	prop_genes_win	prop_genes_genome	enrichment	num_wins	name
9	0.0152801358234295	0.0120013242840589	1.27320414495638	8	locomotion
29	0.0492359932088285	0.0440324449594438	1.11817531945313	17	reproduction
1	0.00169779286926995	0.00173812282734647	0.976796830786644	1	detoxification
21	0.0356536502546689	0.0369144181426916	0.965846193670202	12	immune system process
104	0.176570458404075	0.210561165369972	0.838570864165893	56	response to stimulus
94	0.159592529711375	0.205346796887932	0.777185386526737	52	localization
181	0.307300509337861	0.405975831815925	0.756942865202861	87	cellular process
180	0.305602716468591	0.414666445952657	0.736984435204294	72	metabolic process
109	0.185059422750424	0.257242178447277	0.719397665917191	72	multicellular organismal process
162	0.275042444821732	0.389836119847707	0.705533507077742	81	biological regulation
80	0.135823429541596	0.205098493626883	0.662235139516369	57	developmental process
38	0.0645161290322581	0.107598079788115	0.599602977667494	32	reproductive process
18	0.0305602716468591	0.0533024333719583	0.573337270244335	16	homeostatic process
5	0.00848896434634975	0.0193676543618606	0.438306270224776	5	rhythmic process
1	0.00169779286926995	0.00422115543784142	0.402210459735677	1	biological process involved in interspecies interaction between organisms
1	0.00169779286926995	0.00446945869889091	0.379865434194806	1	pigmentation
2	0.0033955857385399	0.0140705181261381	0.241326275841406	2	growth




select count(distinct gi.ncbi_gene_id) num_genes, 1.0 * count(distinct gi.ncbi_gene_id) / total_genes_windows prop_genes_win, prop_genes_genome, 
       1.0 * count(distinct gi.ncbi_gene_id) / total_genes_windows / prop_genes_genome enrichment, count(distinct dsw_id) num_wins, name
from d_stat_win_50000 w, gene_info gi, sim_go sg, go_edges_tier_2 t, go,
     (select t2.go_id_parent go_id, 1.0 * count(distinct ncbi_gene_id) / total_genes prop_genes_genome
      from sim_go sg2, go_edges_tier_2 t2, (select count(distinct ncbi_gene_id) total_genes from sim_go) x2
      where sg2.go_id = t2.go_id_child
      group by t2.go_id_parent, total_genes) x,
      (select count(distinct ncbi_gene_id) total_genes_windows
       from d_stat_win_50000 w2, gene_info gi2
       where gi2.start between w2.start and w2.end
       and w2.chrom = gi2.chrom
       and d_plus > .104) p
where d_plus > .104
and gi.start between w.start and w.end
and w.chrom = gi.chrom
and gi.ncbi_gene_id = sg.ncbi_gene_id
and sg.go_id = t.go_id_child
and t.go_id_parent = go.go_id
and namespace = 'biological_process'
and x.go_id = go.go_id
group by name, total_genes_windows, prop_genes_genome
order by enrichment desc, num_genes desc, num_wins desc, name;

num_genes	prop_genes_win	prop_genes_genome	enrichment	num_wins	name
1	0.00537634408602151	0.00173812282734647	3.09318996415771	1	detoxification
13	0.0698924731182796	0.0440324449594438	1.5872948500283	7	reproduction
1	0.00537634408602151	0.00422115543784142	1.27366645582964	1	biological process involved in interspecies interaction between organisms
37	0.198924731182796	0.210561165369972	0.944736085750997	19	response to stimulus
6	0.032258064516129	0.0369144181426916	0.873860841892087	4	immune system process
39	0.209677419354839	0.257242178447277	0.815097355419936	21	multicellular organismal process
58	0.311827956989247	0.389836119847707	0.799894984361802	26	biological regulation
59	0.317204301075269	0.414666445952657	0.764962548022235	26	metabolic process
29	0.155913978494624	0.205346796887932	0.759271538964951	17	localization
55	0.295698924731183	0.405975831815925	0.728365832538665	28	cellular process
24	0.129032258064516	0.205098493626883	0.62912338254055	14	developmental process
9	0.0483870967741935	0.107598079788115	0.44970223325062	6	reproductive process
1	0.00537634408602151	0.0120013242840589	0.447979236188357	1	locomotion
4	0.021505376344086	0.0533024333719583	0.40345956054231	4	homeostatic process
1	0.00537634408602151	0.0193676543618606	0.277593971142358	1	rhythmic process


select distinct gi.chrom, gi.start, gi.ncbi_gene_id, synonyms
from d_stat_win_50000 w, gene_info gi, sim_go sg, go_edges_tier_2 t, go, gene_xrefs gx, og2genes og, og2genes og2, genes g
where d_plus > .104
and gi.start between w.start and w.end
and w.chrom = gi.chrom
and gi.ncbi_gene_id = sg.ncbi_gene_id
and sg.go_id = t.go_id_child
and t.go_id_parent = go.go_id
and namespace = 'biological_process'
and name = 'detoxification'
and gi.ncbi_gene_id = gx.external_id
and gx.orthodb_gene_id = og.orthodb_gene_id
and og.orthogroup_id = og2.orthogroup_id
and og.orthodb_gene_id != og2.orthodb_gene_id
and og2.orthodb_gene_id = g.orthodb_gene_id
and orthodb_species_id = '7227_0';

chrom	start	ncbi_gene_id	synonyms
2R	7184270	6733772	Marc


select distinct gi.chrom, gi.start, gi.ncbi_gene_id, synonyms
from d_stat_win_50000 w, gene_info gi, sim_go sg, go_edges_tier_2 t, go, gene_xrefs gx, og2genes og, og2genes og2, genes g
where d_plus > .104
and gi.start between w.start and w.end
and w.chrom = gi.chrom
and gi.ncbi_gene_id = sg.ncbi_gene_id
and sg.go_id = t.go_id_child
and t.go_id_parent = go.go_id
and namespace = 'biological_process'
and name = 'reproduction'
and gi.ncbi_gene_id = gx.external_id
and gx.orthodb_gene_id = og.orthodb_gene_id
and og.orthogroup_id = og2.orthogroup_id
and og.orthodb_gene_id != og2.orthodb_gene_id
and og2.orthodb_gene_id = g.orthodb_gene_id
and orthodb_species_id = '7227_0';

chrom	start	ncbi_gene_id	synonyms
X	41306	6726640	Ggt-1
X	41306	6726640	CG17636-RA;Dmel\\CG17636;EG:23E12.1
3L	4239329	6736763	Hexo2
X	41306	6726640	Dmel\\CG4829
3L	4159344	6736742	Cpr5C
3L	4160415	6736743	Cpr5C
3L	4156175	6736741	Cpr5C
3L	4163767	6736744	Cpr5C
X	41306	6726640	Dmel\\CG1492
X	41306	6726640	CG4829;Dmel\\CG4829
3L	4159344	6736742	Cpr30F
3L	4160415	6736743	Cpr30F
3L	4156175	6736741	Cpr30F
3L	4163767	6736744	Cpr30F
3L	4159344	6736742	Cpr30B
3L	4160415	6736743	Cpr30B
3L	4156175	6736741	Cpr30B
3L	4163767	6736744	Cpr30B
3L	4159344	6736742	Cpr31A
3L	4160415	6736743	Cpr31A
3L	4156175	6736741	Cpr31A
3L	4163767	6736744	Cpr31A
3L	4159344	6736742	CG4588
3L	4160415	6736743	CG4588
3L	4156175	6736741	CG4588
3L	4163767	6736744	CG4588
2R	7160036	6733765	lectin-46Ca
2R	7161877	6733766	lectin-46Ca
2R	7160036	6733765	lectin-46Cb
2R	7161877	6733766	lectin-46Cb
2R	5944914	6733558	CG14747;CG42326;Dmel\\CG42326
2R	7165111	6733767	Dmel\\CG34033
3L	4159344	6736742	Cpr76Ba
3L	4160415	6736743	Cpr76Ba
3L	4156175	6736741	Cpr76Ba
3L	4163767	6736744	Cpr76Ba
3L	4159344	6736742	Cpr64Ab
3L	4163767	6736744	Cpr64Ab
3L	4160415	6736743	Cpr64Ab
3L	4156175	6736741	Cpr64Ab
3L	4159344	6736742	Dmel\\CG13670
3L	4160415	6736743	Dmel\\CG13670
3L	4156175	6736741	Dmel\\CG13670
3L	4163767	6736744	Dmel\\CG13670
3L	4159344	6736742	CG7072
3L	4160415	6736743	CG7072
3L	4156175	6736741	CG7072
3L	4163767	6736744	CG7072
3L	4159344	6736742	Cpr76Bc
3L	4160415	6736743	Cpr76Bc
3L	4156175	6736741	Cpr76Bc
3L	4163767	6736744	Cpr76Bc
3L	17520632	6738509	NUCB1
3L	4159344	6736742	Cpr62Bc
3L	4160415	6736743	Cpr62Bc
3L	4156175	6736741	Cpr62Bc
3L	4163767	6736744	Cpr62Bc
3L	4159344	6736742	Cpr66D
3L	4160415	6736743	Cpr66D
3L	4156175	6736741	Cpr66D
3L	4163767	6736744	Cpr66D
3L	4159344	6736742	Cpr64Ad
3L	4163767	6736744	Cpr64Ad
3L	4160415	6736743	Cpr64Ad
3L	4156175	6736741	Cpr64Ad
3L	4159344	6736742	Cpr62Bb
3L	4160415	6736743	Cpr62Bb
3L	4156175	6736741	Cpr62Bb
3L	4163767	6736744	Cpr62Bb
3L	4156175	6736741	Cpr64Aa
3L	4159344	6736742	Cpr64Aa
3L	4160415	6736743	Cpr64Aa
3L	4163767	6736744	Cpr64Aa
3L	4159344	6736742	Cpr76Bd
3L	4160415	6736743	Cpr76Bd
3L	4156175	6736741	Cpr76Bd
3L	4163767	6736744	Cpr76Bd
3L	4159344	6736742	Cpr76Bb
3L	4160415	6736743	Cpr76Bb
3L	4156175	6736741	Cpr76Bb
3L	4163767	6736744	Cpr76Bb
3L	4160415	6736743	Cpr64Ac
3L	4159344	6736742	Cpr64Ac
3L	4156175	6736741	Cpr64Ac
3L	4163767	6736744	Cpr64Ac
3L	4239329	6736763	Hexo1
3L	4159344	6736742	Cpr66Cb
3L	4160415	6736743	Cpr66Cb
3L	4156175	6736741	Cpr66Cb
3L	4163767	6736744	Cpr66Cb
3L	4159344	6736742	Ccp84Ae
3L	4160415	6736743	Ccp84Ae
3L	4156175	6736741	Ccp84Ae
3L	4163767	6736744	Ccp84Ae
3L	4159344	6736742	Ccp84Ag
3L	4160415	6736743	Ccp84Ag
3L	4156175	6736741	Ccp84Ag
3L	4163767	6736744	Ccp84Ag
3L	4159344	6736742	Cpr92A
3L	4160415	6736743	Cpr92A
3L	4156175	6736741	Cpr92A
3L	4163767	6736744	Cpr92A
3L	4159344	6736742	Ccp84Ac
3L	4160415	6736743	Ccp84Ac
3L	4156175	6736741	Ccp84Ac
3L	4163767	6736744	Ccp84Ac
3L	4159344	6736742	Ccp84Ab
3L	4160415	6736743	Ccp84Ab
3L	4156175	6736741	Ccp84Ab
3L	4163767	6736744	Ccp84Ab
3L	4159344	6736742	Ccp84Ad
3L	4160415	6736743	Ccp84Ad
3L	4156175	6736741	Ccp84Ad
3L	4163767	6736744	Ccp84Ad
3L	4159344	6736742	Edg84A
3L	4160415	6736743	Edg84A
3L	4156175	6736741	Edg84A
3L	4163767	6736744	Edg84A
3L	4159344	6736742	Ccp84Af
3L	4160415	6736743	Ccp84Af
3L	4156175	6736741	Ccp84Af
3L	4163767	6736744	Ccp84Af
3L	4159344	6736742	Ccp84Aa
3L	4160415	6736743	Ccp84Aa
3L	4156175	6736741	Ccp84Aa
3L	4163767	6736744	Ccp84Aa
3R	8956798	6727723	CG5329
3R	8957713	6727724	CG5329
3R	8956798	6727723	Dmel\\CG31419
3R	8957713	6727724	Dmel\\CG31419
