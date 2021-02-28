#include <cstdint>

namespace ECM_opt
{
	const size_t size = 270;
	const size_t min = 18;

	uint64_t B2[size] = {
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5118, 6443, 7996, 9782, 11948, 14587, 17790, 21669, 26344, 31964, 38720, 46742, 56317, 67596, 80914, 96464, 114677, 135758, 160272, 188455, 221100, 258767, 302322, 353256, 412188, 480630, 559641, 651664, 757072, 879676, 1019489, 1181152, 1367239, 1578698, 1821589, 2098581, 2413515, 2774067, 3180251, 3643758, 4163963, 4756255, 5424944, 6179117, 7036618, 7993719, 9083050, 10303871, 11699626, 13266573, 15029514, 17011135, 19258397, 21782074, 24613887, 27784889, 31375952, 35402406, 39904379, 44920965, 50620487, 56972289, 64062544, 71879679, 80680232, 90477466, 101322565, 113633619, 127055385, 142234450, 159042528, 177684676, 198353668, 221261770, 246628178, 274701712, 305727705, 340536586, 379031620, 421598791, 468648424, 520621032, 577997719, 642228136, 712139681, 790222237, 876614632, 971564654, 1076537021, 1193714005, 1320949239, 1463052914, 1619536417, 1788995213, 1978106848, 2185935380, 2414199761, 2664080369, 2944540526, 3246752017, 3577962172, 3947057663, 4345032574, 4788249591, 5273681610, 5805617431, 6387660431, 7024932838, 7709506982, 8469334731, 9316136668, 10225350306, 11211650367, 12322053984, 13497935804, 14818326361, 16256716540, 17825389948, 19490321870, 21351454392, 23378655636, 25591448688, 28065863375, 30695506603, 33540946882, 36734070187, 40119619519, 43892385532, 47991760141, 52364109947, 57308515443, 62589181218, 68169588500, 74393702606, 81185766968, 88532401223, 96700074551, 105561779282, 114704550348, 125207885875, 136269731863, 148552678960, 161874764432, 175959330760, 191658088611, 208619275179, 227001945753, 246293749538, 268009533008, 291093170374, 316343577255, 342911537638, 374152803411, 406180159585, 440799674683, 479114658736, 519606367000, 564427750987, 611746414577, 662653433771, 718985709499, 780430636001, 844846729347, 915849695443, 992897757151, 1075938629040, 1162528454040, 1262327030828, 1363040271103, 1479350526198, 1600619142866, 1731648569219, 1872538245450, 2028706318721, 2186800584642, 2369572709399, 2557772726757, 2770050624917, 2989216577657, 3234882404153, 3487937890070, 3780930586136, 4066575252127, 4393563731309, 4749790972671, 5127636067686, 5537668802637, 5976897329126, 6449434233881, 6957587491919, 7503182263341, 8091105166609, 8695298161223, 9399603596351, 10126828606048, 10901408192240, 11772768037045, 12675488162019, 13632513650759, 14678302339654, 15828215149160, 17063322280989, 18296220497083, 19711533374169, 21238574410702, 22822434354472, 24578387508576, 26458802300846, 28414352334071, 30569601460308, 32882670867690, 35377005129064, 37917857139184, 40771121847905, 43948086166951, 47223004592636, 50608502589479, 54332790134732, 58515525425424, 62667649347385, 67221427706638, 72144996814239, 77654648188045, 83305731135713, 89311665170324, 96012120731097, 102846557546325, 110336088506793, 118140584513188, 126922852025520, 136303235968760, 145910219500832, 156293096086597, 167735757196539, 179379741379102, 191411500856004, 206390068792538, 220870206633117, 236810466514708, 253808917448601, 271115187294142, 290462280420027, 312116841866356, 334263479011988, 356011101351604, 381117851008148, 407913227669492, 437907785791879, 468428481167425, 500883872894296, 536230713310580, 573449613009940 };
	uint64_t B1[size] = {
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 168, 209, 257, 312, 379, 459, 556, 672, 810, 975, 1171, 1402, 1676, 1996, 2371, 2805, 3312, 3894, 4569, 5340, 6230, 7252, 8427, 9794, 11366, 13181, 15264, 17679, 20427, 23610, 27220, 31369, 36126, 41500, 47646, 54622, 62515, 71517, 81606, 93077, 105889, 120423, 136766, 155125, 175928, 199054, 225253, 254578, 287855, 325132, 366898, 413652, 466494, 525606, 591672, 665418, 748619, 841507, 945041, 1060218, 1190199, 1334800, 1495659, 1672272, 1870566, 2090575, 2333982, 2608139, 2907239, 3243276, 3614905, 4025792, 4480004, 4981810, 5535781, 6147001, 6820690, 7574231, 8405416, 9321813, 10331684, 11443882, 12668244, 14035634, 15518735, 17173104, 18994916, 20995513, 23196561, 25650710, 28306239, 31266223, 34516886, 38025437, 41932599, 46215490, 50907544, 56040325, 61768897, 67932563, 74671157, 82166078, 90222095, 99178838, 108964352, 119659780, 131337216, 144089999, 157746409, 172887960, 189717364, 207740495, 227322554, 249169457, 272407401, 298278049, 326476704, 357165948, 389784847, 426048878, 465480846, 508382479, 556114567, 606873125, 661865052, 723108625, 788134643, 860302972, 938639921, 1021961789, 1116122799, 1216416315, 1322554041, 1440311674, 1568215564, 1706618418, 1860012647, 2026315021, 2198192154, 2393759924, 2600175187, 2828752478, 3076340962, 3337938043, 3628009953, 3941411348, 4280485896, 4637289275, 5034385990, 5460174866, 5922299713, 6415329897, 6975216409, 7558069301, 8187101359, 8882286179, 9615367457, 10425651152, 11279253464, 12197521064, 13212278737, 14311539941, 15465050987, 16736769342, 18111036606, 19591016883, 21139261146, 22904556127, 24699832210, 26749186076, 28892689324, 31203667909, 33687469685, 36432726118, 39222818051, 42410145832, 45723085173, 49411807681, 53249599233, 57511603747, 61937084976, 66999694858, 71977125462, 77638563842, 83757377307, 90289876894, 97334816076, 104887243223, 112999415071, 121709441256, 131052453335, 141090874503, 151476946076, 163397153165, 175776619457, 188996846715, 203702458958, 218986123155, 235269548543, 252859946152, 272223422455, 292998798539, 313846516173, 337621227962, 363167855664, 389688105651, 418984978249, 450350504472, 482923880137, 518819645976, 557265986365, 598543951177, 640991359870, 688128169148, 740347235882, 794389497947, 850281158826, 911780278567, 980159017427, 1048437724898, 1123478188653, 1204042703710, 1293442345045, 1385627745403, 1483792832798, 1592582496510, 1704416013399, 1825157381241, 1952534170167, 2094219604216, 2245569956308, 2401255732491, 2568252958864, 2752207540170, 2940793573371, 3139272164850, 3373983259646, 3605234888295, 3859935765019, 4131467895783, 4409689582660, 4718134925394, 5060620319012, 5412602483765, 5760301848054, 6158748553791, 6583484238648, 7054950714200, 7538207932090, 8052377242044, 8605439634612, 9191236809001 };
	uint64_t n[size] = {
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 8, 9, 10, 12, 13, 15, 16, 18, 21, 23, 25, 28, 31, 34, 38, 42, 47, 52, 57, 64, 70, 78, 86, 95, 105, 115, 126, 139, 152, 166, 182, 199, 218, 238, 260, 283, 309, 337, 368, 401, 437, 477, 519, 566, 616, 671, 729, 793, 861, 934, 1013, 1099, 1190, 1289, 1395, 1510, 1632, 1763, 1905, 2058, 2220, 2394, 2582, 2788, 3007, 3242, 3497, 3764, 4058, 4368, 4701, 5060, 5446, 5861, 6307, 6787, 7302, 7844, 8426, 9049, 9718, 10434, 11203, 12009, 12890, 13815, 14803, 15862, 16993, 18177, 19470, 20824, 22269, 23850, 25502, 27267, 29153, 31172, 33267, 35562, 38014, 40569, 43363, 46273, 49378, 52688, 56219, 59981, 64095, 68383, 72831, 77693, 82902, 88253, 94162, 100247, 106724, 113622, 121213, 129040, 137371, 146222, 155319, 165313, 175993, 186929, 198949, 211324, 224485, 238861, 253281, 269015, 286312, 304079, 322871, 342880, 363436, 385253, 409972, 434385, 461276, 488888, 518149, 550228, 583039, 617903, 654832, 695488, 736717, 781047, 827603, 878024, 927114, 982344, 1040839, 1100752, 1166265, 1233308, 1306644, 1384529, 1464191, 1547773, 1639626, 1733843, 1833016, 1938005, 2053656, 2166189, 2295451, 2420949, 2559497, 2705518, 2860118, 3017329, 3197054, 3371231, 3564990, 3759098, 3974408, 4191077, 4431779, 4663406, 4940642, 5211357, 5493990, 5795707, 6111649, 6446001, 6798438, 7170019, 7562249, 7974915, 8432488, 8869900, 9354339, 9868535, 10381524, 10947650, 11550301, 12176327, 12812307, 13481644, 14252356, 14997951, 15779357, 16639058, 17505793, 18419311, 19421691, 20436386, 21503493, 22620717, 23864916, 25106382, 26347460, 27720854, 29231962, 30763059, 32282371, 34042273, 35828542, 37691969, 39545700, 41602089, 43776013, 45943933, 48355774, 50842847, 53511557, 56153800, 58932473, 62011676, 65218889, 68446642, 72038232, 75898577, 79343280, 83442694, 87567435, 91905889, 96723119, 101514497, 106249612, 111511881, 117608213, 123432158, 129545608, 135585798, 142315090, 149395500, 156702217, 164450121 };
}
