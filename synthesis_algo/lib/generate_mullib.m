function [mult_lib, minW, maxW, len_mult_lib] = generate_mullib()

% Zhiyuan Yang   May 26 2017
% this routine generates the special multiplier library 
minW = -141; 
maxW = 205;
mult_lib = [-141	2063	1551;
-140	0	0;
-139	0	0;
-138	0	0;
-137	0	0;
-136	0	0;
-135	0	0;
-134	0	0;
-133	1942	979;
-132	1550	516;
-131	1945	926;
-130	1546	500;
-129	1541	485;
-128	19	55;
-127	1633	732;
-126	1753	770;
-125	0	0;
-124	0	0;
-123	0	0;
-122	0	0;
-121	1944	1213;
-120	1698	793;
-119	0	0;
-118	0	0;
-117	0	0;
-116	1903	1222;
-115	2275	1775;
-114	2024	1211;
-113	1966	1218;
-112	1690	806;
-111	1930	1320;
-110	1934	1367;
-109	2393	1992;
-108	2015	1350;
-107	2068	2192;
-106	1998	1664;
-105	1892	1581;
-104	1714	1080;
-103	2108	1847;
-102	2065	1730;
-101	1909	1560;
-100	1714	1064;
-99	2118	1704;
-98	1714	1049;
-97	1714	1034;
-96	1528	565;
-95	1928	1220;
-94	1904	1231;
-93	2119	1854;
-92	1894	1245;
-91	2194	1913;
-90	1959	1663;
-89	1900	1590;
-88	1653	1035;
-87	2099	1843;
-86	1919	1564;
-85	2056	1608;
-84	1748	1077;
-83	1961	1533;
-82	1799	1035;
-81	1799	1020;
-80	1523	550;
-79	1931	1221;
-78	1950	1240;
-77	1977	1573;
-76	1791	1004;
-75	1951	1530;
-74	1741	1015;
-73	1833	1037;
-72	1518	535;
-71	2064	1285;
-70	1784	976;
-69	1840	993;
-68	1513	520;
-67	1843	940;
-66	1508	504;
-65	1504	489;
-64	19	52;
-63	1675	730;
-62	1616	734;
-61	1994	1221;
-60	1691	753;
-59	2017	1313;
-58	1838	1205;
-57	1869	1202;
-56	1594	764;
-55	1955	1333;
-54	1968	1316;
-53	1979	1656;
-52	1790	1084;
-51	2048	1724;
-50	1790	1069;
-49	1790	1053;
-48	1448	530;
-47	1845	1216;
-46	1819	1230;
-45	2008	1655;
-44	1636	1050;
-43	1881	1556;
-42	1730	1070;
-41	1773	1038;
-40	1527	534;
-39	1968	1228;
-38	1758	1018;
-37	1739	1007;
-36	1522	519;
-35	1795	968;
-34	1505	514;
-33	1500	498;
-32	19	49;
-31	1480	698;
-30	1574	726;
-29	1887	1198;
-28	1588	735;
-27	1935	1299;
-26	1716	1064;
-25	1716	1049;
-24	1425	512;
-23	1868	1222;
-22	1674	1021;
-21	1768	1062;
-20	1488	523;
-19	1736	1010;
-18	1456	580;
-17	1511	486;
-16	19	46;
-15	1633	705;
-14	1507	710;
-13	1679	1056;
-12	1462	829;
-11	1636	1013;
-10	1430	511;
-9	1470	504;
-8	19	43;
-7	1583	699;
-6	1506	761;
-5	1366	527;
-4	19	40;
-3	1302	679;
-2	19	37;
-1	19	35;
0	0	3;
1	0	26;
2	0	29;
3	1289	704;
4	0	32;
5	1344	505;
6	1361	766;
7	1529	683;
8	0	35;
9	1470	496;
10	1430	503;
11	1636	1004;
12	1490	830;
13	1679	1047;
14	1557	706;
15	1633	697;
16	0	38;
17	1511	478;
18	1456	571;
19	1736	1002;
20	1488	514;
21	1768	1053;
22	1674	1012;
23	1868	1214;
24	1425	504;
25	1716	1040;
26	1716	1056;
27	1935	1290;
28	1588	727;
29	1887	1190;
30	1574	718;
31	1480	689;
32	0	41;
33	1500	490;
34	1505	505;
35	1795	960;
36	1522	511;
37	1739	999;
38	1758	1009;
39	1968	1220;
40	1527	526;
41	1773	1030;
42	1730	1062;
43	1881	1547;
44	1636	1041;
45	2008	1647;
46	1819	1221;
47	1845	1207;
48	1448	522;
49	1790	1045;
50	1790	1060;
51	2048	1716;
52	1790	1075;
53	1979	1648;
54	1968	1307;
55	1955	1325;
56	1594	756;
57	1869	1194;
58	1838	1197;
59	2017	1305;
60	1691	745;
61	1994	1212;
62	1616	726;
63	1675	721;
64	0	44;
65	1504	481;
66	1508	496;
67	1843	801;
68	1513	511;
69	1840	985;
70	1784	968;
71	2064	1276;
72	1518	526;
73	1833	1028;
74	1741	1007;
75	1951	1521;
76	1791	996;
77	1977	1565;
78	1950	1231;
79	1931	1213;
80	1523	542;
81	1799	1011;
82	1799	1026;
83	1961	1525;
84	1748	1069;
85	2056	1600;
86	1919	1556;
87	2099	1835;
88	1653	1027;
89	1900	1581;
90	1959	1654;
91	2194	1905;
92	1894	1236;
93	2119	1846;
94	1904	1222;
95	1928	1211;
96	1528	557;
97	1714	1025;
98	1714	1041;
99	2118	1696;
100	1714	1056;
101	1909	1552;
102	2065	1721;
103	2108	1839;
104	1714	1071;
105	1892	1573;
106	1998	1655;
107	2068	2183;
108	2015	1342;
109	2393	1984;
110	1934	1358;
111	1930	1311;
112	1690	797;
113	0	0;
114	2024	1203;
115	2275	1766;
116	0	0;
117	0	0;
118	0	0;
119	2085	1310;
120	0	0;
121	0	0;
122	0	0;
123	0	0;
124	0	0;
125	2018	1234;
126	1753	762;
127	1633	723;
128	0	0;
129	0	0;
130	1546	492;
131	0	0;
132	0	0;
133	1942	971;
134	0	0;
135	0	0;
136	1555	522;
137	0	0;
138	0	0;
139	2051	1499;
140	0	0;
141	0	0;
142	0	0;
143	0	0;
144	0	0;
145	1824	997;
146	0	0;
147	0	0;
148	0	0;
149	0	0;
150	0	0;
151	0	0;
152	0	0;
153	0	0;
154	2016	1573;
155	0	0;
156	0	0;
157	0	0;
158	0	0;
159	0	0;
160	0	0;
161	0	0;
162	0	0;
163	0	0;
164	0	0;
165	0	0;
166	0	0;
167	0	0;
168	0	0;
169	0	0;
170	0	0;
171	0	0;
172	0	0;
173	0	0;
174	0	0;
175	0	0;
176	0	0;
177	0	0;
178	0	0;
179	0	0;
180	0	0;
181	0	0;
182	0	0;
183	0	0;
184	0	0;
185	0	0;
186	0	0;
187	0	0;
188	0	0;
189	0	0;
190	0	0;
191	0	0;
192	0	0;
193	0	0;
194	0	0;
195	0	0;
196	0	0;
197	0	0;
198	0	0;
199	0	0;
200	0	0;
201	0	0;
202	0	0;
203	0	0;
204	0	0;
205	2114	2059];
        
  len_mult_lib = [0     0 0;      % dummy: for non len_mult case
                  52	121	0.002125907;
1604	863	0.023642921;
1738	1559	0.06236947;
1986	2266	0.11211913;
2604	3170	0.233218952;
2145	3675	0.228095959;
2947	4638	0.383347236;
3070	5750	0.496439068;
3215	6050	0.550806076;
3298	7237	0.652021993;
3399	7375	0.694912079;
3666	8690	0.844003428;
3737	8824	0.895502757;
3883	10013	1.010393074;
3989	10051	1.056710871
                  ];
              
   power = [0.117350561
0.0667826
0.119699685
0.068643411
0.067675266
0.029205409
0.076551228
0.064577001
0.064478875
0.028400176
0.059963718
0.027469265
0.022994908
0.000104117
0.031825933
0.033472093
0.07423258
0.036961782
0.080573814
0.077394737
0.077948077
0.039768105
0.083243267
0.083503603
0.139930398
0.079125213
0.130749014
0.080933584
0.081119656
0.039076216
0.082403448
0.08673908
0.158793368
0.085618373
0.178104511
0.137621583
0.131482306
0.0721491
0.143540719
0.130867315
0.130725678
0.072500868
0.127330318
0.071788137
0.070099843
0.02600194
0.07737433
0.077431101
0.142005191
0.078358566
0.150525423
0.125831651
0.122940963
0.070633941
0.137627755
0.126144853
0.137185692
0.072362085
0.122360502
0.071245242
0.07297832
0.027456373
0.078218241
0.079412705
0.119537448
0.06868033
0.121326545
0.07081585
0.072744728
0.027922918
0.080083365
0.065297203
0.067002038
0.028410879
0.062880173
0.027599116
0.026708079
0.000104117
0.03184132
0.034204194
0.075241573
0.038351111
0.08070691
0.075985671
0.077599759
0.037450812
0.083652075
0.082865721
0.135216258
0.070829312
0.128765367
0.071098607
0.070390915
0.024499143
0.074860634
0.076100119
0.123327353
0.069268713
0.122628401
0.070513873
0.069669637
0.026141977
0.077103688
0.067544878
0.0688063
0.026525214
0.063394497
0.027000355
0.026430522
0.000104117
0.03192449
0.035834118
0.074169118
0.035676333
0.080370246
0.068365895
0.068148732
0.022980302
0.074225184
0.067426418
0.068406318
0.024960879
0.065570879
0.026931223
0.025003073
0.000104117
0.034287212
0.033468616
0.065735011
0.033888112
0.06473886
0.023480857
0.023789963
0.000104117
0.031607539
0.031702397
0.022407145
0.000104117
0.027304882
0.000104117
0.000104117
0
0
0
0.028533419
0
0.022095783
0.031447127
0.031404331
0
0.023685846
0.023308568
0.064634743
0.033630602
0.065630894
0.033134475
0.034183095
0
0.024898956
0.026827106
0.065466762
0.024856762
0.068302201
0.067322301
0.074121067
0.022876185
0.068044615
0.068261778
0.080266129
0.035572216
0.074065001
0.035730001
0.031820373
0
0.026326405
0.026896238
0.06329038
0.026421097
0.068702183
0.067440761
0.076999571
0.02603786
0.06956552
0.070409756
0.122524284
0.069164596
0.123223236
0.075996002
0.074756517
0.024395026
0.070286798
0.07099449
0.12866125
0.070725195
0.135112141
0.082761604
0.083547958
0.037346695
0.077495642
0.075881554
0.080602793
0.038246994
0.075137456
0.034100077
0.031737203
0.000104117
0.026708079
0.027599116
0.062880173
0.028410879
0.067002038
0.065297203
0.080083365
0.027922918
0.072744728
0.070711733
0.121222428
0.068576213
0.119433331
0.079308588
0.078114124
0.027456373
0.07297832
0.071245242
0.122360502
0.072362085
0.137185692
0.126144853
0.137627755
0.070633941
0.122940963
0.125727534
0.150421306
0.078254449
0.141901074
0.077326984
0.077270213
0.02600194
0.070099843
0.071788137
0.127330318
0.072500868
0.130725678
0.130867315
0.143540719
0.0721491
0.131482306
0.137517466
0.178000394
0.085514256
0.158689251
0.086634963
0.082299331
0.039076216
0
0.080829467
0.130644897
0
0
0
0.08313915
0
0
0
0
0
0.074128463
0.033367976
0.031721816
0
0
0.027365148
0
0
0.064374758
0
0.076447111
0.029101292
0
0
0.119595568
0
0
0
0
0
0.070321625
0
0
0
0
0
0
0
0
0.122446254
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0.164892811
];
mult_lib = [mult_lib , power];              
              