import networkx as nx
import matplotlib.pyplot as plt
import numpy as np

L = [(1,2),
(2,3),
(3,4),
(4,5),
(5,4),
(4,6),
(6,4),
(4,7),
(7,8),
(8,9),
(9,10),
(10,11),
(11,12),
(12,13),
(13,14),
(14,1),
(1,4),
(4,15),
(15,16),
(16,17),
(17,18),
(18,19),
(19,20),
(20,4),
(4,1),
(8,21),
(21,22),
(22,23),
(23,20),
(23,24),
(24,1),
(23,25),
(25,13),
(22,26),
(26,12),
(7,27),
(27,28),
(28,10),
(1,29),
(29,30),
(30,29),
(29,28),
(6,31),
(31,30),
(19,32),
(32,33),
(33,14),
(4,31),
(8,34),
(34,12)]
tuplesmk1 = L
Gmk1 = nx.DiGraph()
xmk1 = list(range(1,34))
Gmk1.add_nodes_from(xmk1)
Gmk1.add_edges_from(tuplesmk1)
#plt.plot()
#nx.draw(Gmk1, with_labels=True, font_weight='bold')
#plt.show()
Amk1 = nx.to_numpy_matrix(Gmk1)
smk1 = Amk1[0].sum()
Smk1list = [smk1];
Atot = Amk1
for i in range(15):
   Aux = np.dot(Amk1,Atot)
   saux = Aux[0].sum()
   smk1 = smk1 + saux
   Smk1list .append(smk1)
   Atot = Aux
#print('Linguagem k= 1 monolitico')
#print(Smk1list)
#print(' ')

L2 = [
(1,	2),
(2,	3),
(3,	4),
(4,	5),
(5,	6),
(6,	7),
(7,	8),
(8,	9),
(9,	10),
(10,	11),
(11,	12),
(12,	13),
(13,	14),
(14,	15),
(15,	16),
(16,	17),
(17,	18),
(18,	5),
(6,	19),
(19,	20),
(20,	21),
(21,	22),
(22,	23),
(23,	24),
(24,	25),
(25,	26),
(26,	18),
(26,	2),
(10,	27),
(27,	28),
(28,	29),
(29,	30),
(30,	25),
(29,	31),
(31,	32),
(32,	2),
(32,	18),
(29,	33),
(33,	34),
(34,	16),
(17,	2),
(28,	35),
(35,	36),
(36,	15),
(9,	37),
(37,	38),
(38,	39),
(39,	13),
(6,	26),
(26,	40),
(40,	41),
(41,	42),
(42,	43),
(43,	39),
(7,	44),
(44,	45),
(45,	42),
(23,	46),
(46,	47),
(47,	48),
(48,	17),
(6,	49),
(49,	45),
(10,	50),
(50,	51),
(51,	15)]

tuplesmk2 = L2
Gmk2 = nx.DiGraph()
xmk2 = list(range(1,51))
Gmk2.add_nodes_from(xmk2)
Gmk2.add_edges_from(tuplesmk2)
#plt.plot()
#nx.draw(Gmk1, with_labels=True, font_weight='bold')
#plt.show()
Amk2 = nx.to_numpy_matrix(Gmk2)
smk2 = Amk2[0].sum()
Smk2list = [smk2];
Atot = Amk2
for i in range(15):
   Aux = np.dot(Amk2,Atot)
   saux = Aux[0].sum()
   smk2 = smk2 + saux
   Smk2list .append(smk2)
   Atot = Aux
   
 



L3 = [
(1,	2),
(2,	3),
(3,	4),
(4,	5),
(5,	6),
(6,	7),
(7,	8),
(8,	9),
(9,	10),
(10,	11),
(11,	12),
(12,	13),
(13,	14),
(14,	15),
(15,	16),
(16,	17),
(17,	18),
(18,	19),
(19,	6),
(6,	20),
(20,	21),
(21,	22),
(22,	23),
(23,	24),
(24,	25),
(25,	26),
(26,	27),
(27,	28),
(28,	19),
(27,	29),
(29,	3),
(10,	30),
(30,	31),
(31,	32),
(32,	33),
(33,	34),
(34,	27),
(32,	35),
(35,	36),
(36,	37),
(37,	3),
(36,	38),
(38,	19),
(32,	39),
(39,	40),
(40,	41),
(41,	17),
(17,	42),
(42,	3),
(31,	43),
(43,	44),
(44,	45),
(45,	16),
(9,	46),
(46,	47),
(47,	48),
(48,	49),
(49,	14),
(6,	50),
(50,	51),
(51,	52),
(52,	53),
(53,	54),
(54,	55),
(55,	49),
(7,	56),
(56,	57),
(57,	58),
(58,	54),
(24,	59),
(59,	60),
(60,	61),
(61,	62),
(62,	18),
(6,	63),
(63,	64),
(64,	58),
(10,	65),
(65,	66),
(66,	67),
(67,	16)]

tuplesmk3 = L3
Gmk3 = nx.DiGraph()
xmk3 = list(range(1,67))
Gmk3.add_nodes_from(xmk3)
Gmk3.add_edges_from(tuplesmk3)
#plt.plot()
#nx.draw(Gmk1, with_labels=True, font_weight='bold')
#plt.show()
Amk3 = nx.to_numpy_matrix(Gmk3)
smk3 = Amk3[0].sum()
Smk3list = [smk3];
Atot = Amk3
for i in range(15):
   Aux = np.dot(Amk3,Atot)
   saux = Aux[0].sum()
   smk3 = smk3 + saux
   Smk3list .append(smk3)
   Atot = Aux
print('Linguagem k= 3 monolitico')
print(Smk3list)
print(' ')     
 

L4 = [(	1	,	2	)	,
(	2	,	3	)	,
(	3	,	4	)	,
(	4	,	5	)	,
(	5	,	6	)	,
(	6	,	7	)	,
(	7	,	8	)	,
(	8	,	9	)	,
(	9	,	10	)	,
(	10	,	11	)	,
(	11	,	12	)	,
(	12	,	13	)	,
(	13	,	14	)	,
(	14	,	15	)	,
(	15	,	16	)	,
(	16	,	17	)	,
(	17	,	18	)	,
(	18	,	19	)	,
(	19	,	20	)	,
(	20	,	21	)	,
(	21	,	22	)	,
(	22	,	23	)	,
(	23	,	24	)	,
(	24	,	25	)	,
(	25	,	26	)	,
(	26	,	27	)	,
(	27	,	28	)	,
(	28	,	29	)	,
(	29	,	30	)	,
(	30	,	20	)	,
(	28	,	31	)	,
(	31	,	32	)	,
(	32	,	4	)	,
(	10	,	33	)	,
(	33	,	34	)	,
(	34	,	35	)	,
(	35	,	36	)	,
(	36	,	37	)	,
(	37	,	38	)	,
(	38	,	31	)	,
(	35	,	39	)	,
(	39	,	40	)	,
(	40	,	41	)	,
(	41	,	42	)	,
(	42	,	4	)	,
(	40	,	43	)	,
(	43	,	44	)	,
(	44	,	20	)	,
(	35	,	45	)	,
(	45	,	46	)	,
(	46	,	47	)	,
(	47	,	48	)	,
(	48	,	18	)	,
(	38	,	29	)	,
(	48	,	49	)	,
(	49	,	50	)	,
(	50	,	4	)	,
(	34	,	51	)	,
(	51	,	52	)	,
(	52	,	53	)	,
(	53	,	54	)	,
(	54	,	17	)	,
(	17	,	49	)	,
(	9	,	55	)	,
(	55	,	56	)	,
(	56	,	57	)	,
(	57	,	58	)	,
(	58	,	59	)	,
(	59	,	15	)	,
(	6	,	60	)	,
(	60	,	61	)	,
(	61	,	62	)	,
(	62	,	63	)	,
(	63	,	64	)	,
(	64	,	65	)	,
(	65	,	66	)	,
(	66	,	59	)	,
(	7	,	67	)	,
(	67	,	68	)	,
(	68	,	69	)	,
(	69	,	70	)	,
(	70	,	65	)	,
(	25	,	71	)	,
(	71	,	72	)	,
(	72	,	73	)	,
(	73	,	74	)	,
(	74	,	75	)	,
(	75	,	19	)	,
(	6	,	76	)	,
(	76	,	77	)	,
(	77	,	78	)	,
(	78	,	70	)	,
(	10	,	79	)	,
(	79	,	80	)	,
(	80	,	81	)	,
(	81	,	82	)	,
(	82	,	17	) ]	

tuplesmk4 = L4
Gmk4 = nx.DiGraph()
xmk4 = list(range(1,82))
Gmk4.add_nodes_from(xmk4)
Gmk4.add_edges_from(tuplesmk4)
Amk4 = nx.to_numpy_matrix(Gmk4)
smk4 = Amk4[0].sum()
Smk4list = [smk4];
Atot = Amk4
for i in range(15):
   Aux = np.dot(Amk4,Atot)
   saux = Aux[0].sum()
   smk4 = smk4 + saux
   Smk4list .append(smk4)
   Atot = Aux
print('Linguagem k= 4 monolitico')
print(Smk4list)
print(' ') 







LP2 = [
(1,	2),
(52,	2),
(59,	2),
(63,	2),
(66,	2),
(2,	3),
(3,	4),
(4,	5),
(57,	5),
(5,	6),
(6,	7),
(6,	8),
(6,	9),
(30,	9),
(6,	10),
(30,	10),
(6,	11),
(20,	11),
(53,	11),
(7,	12),
(21,	12),
(28,	12),
(7,	13),
(21,	13),
(28,	13),
(7,	14),
(21,	14),
(28,	14),
(7,	15),
(21,	15),
(28,	15),
(7,	16),
(21,	16),
(28,	16),
(8,	17),
(23,	17),
(29,	17),
(8,	18),
(23,	18),
(29,	18),
(9,	19),
(9,	20),
(15,	20),
(24,	20),
(9,	21),
(10,	22),
(10,	23),
(11,	24),
(11,	25),
(11,	26),
(11,	27),
(11,	28),
(11,	29),
(11,	30),
(14,	30),
(19,	30),
(52,	30),
(59,	30),
(63,	30),
(66,	30),
(11,	31),
(14,	31),
(19,	31),
(52,	31),
(59,	31),
(12,	32),
(20,	32),
(13,	33),
(37,	33),
(14,	34),
(19,	34),
(14,	35),
(19,	35),
(14,	36),
(19,	36),
(16,	37),
(25,	37),
(40,	37),
(20,	38),
(22,	39),
(31,	40),
(31,	41),
(33,	42),
(35,	42),
(44,	42),
(36,	43),
(38,	43),
(36,	44),
(38,	44),
(39,	45),
(42,	46),
(45,	47),
(46,	48),
(47,	49),
(47,	50),
(48,	51),
(49,	52),
(49,	53),
(60,	53),
(49,	54),
(50,	55),
(51,	56),
(54,	57),
(65,	57),
(54,	58),
(54,	59),
(65,	59),
(54,	60),
(65,	60),
(54,	61),
(55,	62),
(56,	63),
(62,	63),
(56,	64),
(61,	64),
(62,	64),
(56,	65),
(62,	65),
(64,	66)]

tuplespk2 = LP2
Gpk2 = nx.DiGraph()
xpk2 = list(range(1,66))
Gpk2.add_nodes_from(xpk2)
Gpk2.add_edges_from(tuplespk2)
#plt.plot()
#nx.draw(Gmk1, with_labels=True, font_weight='bold')
#plt.show()
Apk2 = nx.to_numpy_matrix(Gpk2)
spk2 = Apk2[0].sum()
Spk2list = [spk2];
Atot = Apk2
for i in range(15):
   Aux = np.dot(Apk2,Atot)
   saux = Aux[0].sum()
   spk2 = spk2 + saux
   Spk2list.append(spk2)
   Atot = Aux

   
print('Linguagem k= 2 paralelo')
print(Spk2list)
print(' ')


print('Linguagem k= 2 monolitico')
print(Smk2list)
print(' ')

Lexc2 = np.subtract(Spk2list, Smk2list)

print('Linguagem k= 2 em excesso')
print(Lexc2)
print(' ')

LP3 = [
(1,	2),
(2,	3),
(2,	3),
(117,	3),
(136,	3),
(139,	3),
(145,	3),
(152,	3),
(2,	4),
(3,	5),
(6,	5),
(4,	6),
(5,	7),
(7,	8),
(135,	8),
(138,	8),
(151,	8),
(8,	9),
(8,	10),
(8,	11),
(8,	12),
(8,	13),
(9,	14),
(9,	15),
(9,	16),
(10,	17),
(11,	18),
(11,	19),
(11,	20),
(12,	21),
(12,	22),
(13,	23),
(13,	24),
(13,	25),
(13,	26),
(13,	27),
(38,	27),
(14,	28),
(39,	28),
(45,	28),
(40,	29),
(46,	29),
(15,	30),
(40,	30),
(46,	30),
(15,	31),
(40,	31),
(46,	31),
(16,	32),
(41,	32),
(47,	32),
(18,	33),
(18,	34),
(18,	35),
(19,	36),
(19,	37),
(19,	38),
(20,	39),
(20,	40),
(20,	41),
(21,	42),
(22,	43),
(23,	44),
(25,	45),
(25,	46),
(25,	47),
(26,	48),
(27,	49),
(27,	50),
(28,	51),
(54,	51),
(55,	51),
(29,	52),
(33,	52),
(59,	52),
(30,	53),
(34,	53),
(60,	53),
(31,	54),
(35,	54),
(32,	55),
(44,	55),
(62,	55),
(36,	56),
(37,	57),
(37,	58),
(38,	59),
(38,	60),
(42,	61),
(49,	62),
(51,	63),
(52,	63),
(64,	63),
(53,	64),
(56,	64),
(68,	64),
(57,	65),
(57,	66),
(57,	67),
(58,	68),
(61,	69),
(63,	70),
(101,	70),
(65,	71),
(65,	72),
(65,	73),
(66,	74),
(66,	75),
(66,	76),
(67,	77),
(67,	78),
(67,	79),
(69,	80),
(69,	81),
(70,	82),
(87,	82),
(97,	82),
(71,	83),
(94,	83),
(98,	83),
(72,	84),
(95,	84),
(99,	84),
(72,	85),
(95,	85),
(99,	85),
(72,	86),
(95,	86),
(99,	86),
(73,	87),
(96,	87),
(100,	87),
(74,	88),
(74,	89),
(74,	90),
(75,	91),
(75,	92),
(75,	93),
(76,	94),
(76,	95),
(76,	96),
(77,	97),
(78,	98),
(78,	99),
(78,	100),
(79,	101),
(80,	102),
(80,	103),
(80,	104),
(81,	105),
(82,	106),
(83,	106),
(109,	106),
(84,	107),
(88,	107),
(114,	107),
(85,	108),
(89,	108),
(115,	108),
(86,	109),
(90,	109),
(116,	109),
(91,	110),
(92,	111),
(92,	112),
(92,	113),
(93,	114),
(93,	115),
(93,	116),
(102,	117),
(128,	117),
(102,	118),
(128,	118),
(103,	119),
(129,	119),
(104,	120),
(130,	120),
(104,	121),
(130,	121),
(104,	122),
(130,	122),
(105,	123),
(106,	124),
(107,	124),
(127,	124),
(106,	125),
(107,	125),
(127,	125),
(106,	126),
(107,	126),
(127,	126),
(108,	127),
(110,	127),
(134,	127),
(111,	128),
(111,	129),
(111,	130),
(112,	131),
(112,	132),
(112,	133),
(113,	134),
(118,	135),
(137,	135),
(140,	135),
(146,	135),
(153,	135),
(119,	136),
(119,	137),
(120,	138),
(131,	138),
(148,	138),
(121,	139),
(132,	139),
(149,	139),
(121,	140),
(132,	140),
(149,	140),
(122,	141),
(133,	141),
(150,	141),
(123,	142),
(123,	143),
(123,	144),
(124,	145),
(124,	146),
(142,	146),
(125,	147),
(143,	147),
(126,	148),
(144,	148),
(126,	149),
(144,	149),
(126,	150),
(144,	150),
(141,	151),
(147,	152),
(147,	153)]

tuplespk3 = LP3
Gpk3 = nx.DiGraph()
xpk3 = list(range(1,153))
Gpk3.add_nodes_from(xpk3)
Gpk3.add_edges_from(tuplespk3)
#plt.plot()
#nx.draw(Gmk1, with_labels=True, font_weight='bold')
#plt.show()
Apk3 = nx.to_numpy_matrix(Gpk3)
spk3 = Apk3[0].sum()
Spk3list = [spk3];
Atot = Apk3
for i in range(20):
   Aux = np.dot(Apk3,Atot)
   saux = Aux[0].sum()
   spk3 = spk3 + saux
   Spk3list.append(spk3)
   Atot = Aux
print(Spk3list)




