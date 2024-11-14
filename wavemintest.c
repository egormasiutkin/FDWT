#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include  "wavemin.h"


int main() {
	wave_object obj;
	wt_object wt;
	float *inp, *out;
	int N, i, J;
	float temp[2560] = { -59.7950, 4.0439, 33.9908, 30.3984, 20.0055, 21.0011, 28.1391, 30.3272, 28.0846, 27.2526, 26.0691, 19.4085, 11.4115, 12.5601, 23.5836, 31.9126, 28.1999, 19.3337, 18.9559, 29.2465, 39.1668, 39.6193, 33.5647, 29.5197, 30.4187, 32.4960, 32.1013, 29.4412, 26.3245, 23.5353, 20.9678, 18.3656, 15.3812, 12.1440, 9.7253, 8.5030, 6.2194, 0.1039, -8.0765, -11.5622, -5.7367, 5.3937, 12.8365, 13.0245, 11.8923, 17.0998, 28.1716, 36.7837, 36.2417, 28.2826, 20.0223, 16.2011, 15.8355, 15.2983, 12.6967, 9.0584, 6.6505, 6.9952, 10.0024, 14.2801, 18.1956, 20.7979, 21.7371, 20.6550, 17.4827, 13.5954, 11.6474, 13.0241, 15.5870, 15.2115, 10.1893, 3.4121, -0.4752, 0.4157, 4.2209, 8.1496, 10.9058, 12.3083, 12.3053, 11.6943, 12.8417, 17.7006, 24.7652, 29.5582, 29.2801, 26.2541, 25.5880, 29.6682, 36.0756, 40.9545, 43.1826, 44.6349, 47.4262, 52.1109, 58.1899, 65.0639, 72.2348, 79.5909, 88.0203, 98.9328, 112.3937, 126.2235, 137.7784, 146.2848, 152.4715, 156.2937, 156.8363, 155.3049, 156.4006, 164.2337, 176.5193, 184.8836, 182.8141, 173.0712, 165.7518, 168.3136, 178.0206, 184.7646, 181.0982, 169.0514, 157.2349, 151.8717, 151.0426, 147.5016, 136.5964, 121.0394, 107.9118, 101.3429, 98.5368, 93.0557, 81.6481, 67.4840, 56.7407, 52.3865, 51.4792, 48.3141, 39.9587, 28.8311, 20.1992, 17.2160, 18.0757, 17.6893, 12.1384, 1.7789, -9.6330, -18.5549, -24.3995, -29.1870, -35.1893, -43.1334, -52.0986, -60.3458, -66.1886, -68.8239, -68.9063, -68.0416, -67.1862, -65.6841, -62.4716, -58.2039, -55.3053, -55.2648, -56.3089, -54.6454, -48.2687, -39.2018, -31.9375, -29.9988, -33.8297, -41.0728, -48.2735, -52.9336, -54.9686, -56.2811, -58.2447, -59.5842, -57.6633, -52.2918, -47.2487, -46.8044, -50.9582, -55.3983, -56.5941, -55.8984, -57.6965, -63.9831, -71.8104, -76.1826, -74.7822, -69.7838, -65.7999, -66.2694, -70.7694, -74.7539, -72.3975, -61.3237, -45.4103, -32.2463, -26.7320, -26.9486, -26.6561, -21.7809, -14.0725, -8.6248, -8.6898, -13.2697, -18.9655, -23.1840, -25.6544, -27.7698, -30.8531, -34.5347, -36.4765, -34.2818, -28.3294, -22.3262, -20.1608, -22.2108, -25.5425, -27.6640, -28.6569, -28.9050, -26.5169, -19.6634, -11.5258, -9.8165, -18.6748, -31.9704, -37.8867, -31.5748, -21.4192, -20.5189, -32.3447, -46.9815, -52.2261, -46.5895, -39.5891, -40.0247, -46.4237, -50.0476, -45.6464, -37.0391, -32.0683, -34.1185, -40.0730, -45.5832, -49.9898, -55.1939, -61.3358, -65.7847, -66.6535, -65.6661, -66.2505, -69.3016, -72.0770, -71.4378, -67.2015, -61.7666, -57.1379, -53.1753, -48.5257, -42.3704, -35.0608, -27.8033, -22.2380, -19.7822, -20.4801, -22.3834, -22.7432, -20.3440, -16.5514, -13.7662, -12.8835, -12.5007, -10.6333, -6.9933, -3.4251, -2.1108, -3.3848, -5.3387, -5.6615, -3.9491, -2.1165, -2.3498, -4.6694, -6.8142, -6.7302, -4.9186, -3.8467, -4.9520, -6.7551, -6.5169, -3.7083, -1.1584, -2.2015, -6.8105, -11.1058, -10.9556, -5.8284, 1.0972, 6.2828, 8.7552, 9.8080, 10.6720, 11.3068, 11.4272, 11.8168, 13.6784, 16.4211, 16.6554, 10.4583, -2.4402, -16.7012, -24.6881, -23.0940, -16.6456, -14.3582, -21.1504, -33.1150, -41.5197, -41.3794, -35.5737, -30.7501, -30.3454, -32.3146, -32.7911, -30.4567, -27.0974, -25.1056, -25.3322, -26.7811, -27.3014, -24.5889, -17.6108, -7.8754, 0.9127, 5.1861, 4.0949, 0.1621, -2.4653, -0.8271, 4.8403, 11.2522, 14.5206, 13.2327, 9.4851, 6.7720, 6.8129, 8.5064, 9.6206, 8.8154, 5.9763, 1.5657, -3.4610, -7.5502, -9.6525, -10.8539, -13.8787, -19.7795, -25.2560, -24.5960, -15.4452, -2.6715, 4.4332, -0.1073, -13.3332, -25.4255, -27.8956, -20.4104, -10.0725, -4.3509, -4.8426, -7.4973, -7.8213, -4.9472, -1.0135, 1.9169, 3.3379, 3.4512, 2.3614, 0.5889, -0.5536, -0.0526, 1.5205, 2.5851, 2.5492, 2.5811, 3.9142, 5.7282, 5.4431, 1.4816, -4.4181, -7.9640, -5.6674, 2.2480, 12.0091, 19.3439, 22.6814, 23.7164, 25.0432, 27.3346, 28.8690, 27.6493, 23.5506, 18.1013, 12.5721, 7.2175, 2.7441, 1.4652, 5.3653, 12.6854, 18.0026, 17.4441, 13.6370, 13.6802, 21.2198, 31.6505, 36.7119, 33.5417, 27.5146, 25.6314, 28.6586, 31.2351, 28.9788, 23.5847, 20.4869, 22.8925, 29.0685, 34.6100, 36.1358, 33.1489, 27.8103, 23.4829, 22.5778, 24.8737, 28.0199, 30.0212, 30.9420, 31.8737, 32.8797, 32.9509, 31.8830, 31.0326, 31.6126, 32.9860, 33.4744, 32.4164, 30.3713, 27.3243, 21.9797, 13.6991, 4.4604, -2.3466, -5.2704, -6.1055, -7.3566, -8.9643, -8.3861, -3.9881, 2.3166, 6.3848, 5.8908, 2.9720, 2.7499, 8.9137, 20.1928, 30.8302, 34.8559, 30.7338, 22.5139, 16.5504, 16.6611, 21.7722, 27.5797, 30.2312, 28.8976, 25.7791, 24.2349, 26.4876, 32.1441, 38.4919, 42.7792, 44.9788, 48.3045, 56.2123, 68.2606, 79.1964, 83.0144, 78.5763, 71.1958, 68.2528, 72.9028, 81.9509, 89.6737, 93.3143, 94.6340, 96.5622, 99.3084, 100.3702, 97.9593, 93.4044, 89.6031, 87.4778, 84.7557, 78.6122, 69.2496, 60.4286, 56.3010, 57.8119, 62.2976, 66.4799, 69.5246, 72.9154, 77.5842, 81.9893, 83.7742, 83.1662, 83.6577, 88.3838, 96.0082, 101.2348, 100.0918, 94.4092, 90.3465, 92.2631, 98.3850, 102.7629, 101.2355, 95.0863, 88.9674, 85.9189, 85.0095, 83.6679, 81.3932, 80.3236, 82.2012, 85.6320, 87.0321, 84.0922, 77.6718, 70.2262, 63.1564, 56.3460, 49.8298, 44.6951, 41.7037, 39.6032, 35.8159, 29.1402, 21.3326, 15.5729, 13.4102, 13.4555, 13.0240, 10.7487, 7.4817, 4.9688, 4.1120, 4.6265, 5.9656, 7.9857, 10.5246, 12.7553, 13.5662, 12.7338, 11.3350, 10.6512, 10.8535, 11.1827, 11.3189, 11.8699, 13.0880, 13.6909, 11.9960, 8.4685, 6.2643, 8.2991, 13.8353, 18.7151, 19.2671, 15.5518, 10.5180, 6.5943, 3.8713, 1.3746, -1.0208, -2.5739, -2.9648, -2.7200, -2.4438, -2.4810, -3.6490, -7.5162, -14.9884, -24.5701, -32.7663, -36.7826, -36.8179, -35.5215, -35.2771, -36.1835, -36.6004, -35.3716, -33.3259, -32.6006, -34.5462, -38.3410, -41.7426, -43.2193, -43.1817, -42.9402, -42.5829, -40.3239, -34.2755, -24.7980, -14.8332, -7.7484, -4.7331, -4.0685, -2.7169, 1.2865, 7.6105, 14.4403, 20.4127, 25.6919, 31.0201, 35.8724, 38.1721, 36.3573, 31.5391, 26.9547, 25.0496, 25.5124, 26.4220, 26.7577, 27.0700, 27.9251, 28.7295, 28.5855, 27.6339, 26.7146, 25.8847, 24.1166, 20.6907, 16.2595, 12.1340, 8.9621, 6.5814, 4.7741, 3.4643, 2.3297, 1.0468, 0.0784, 0.3783, 1.8615, 3.0097, 2.9534, 3.4253, 7.1971, 14.2642, 20.8478, 23.0997, 21.0160, 17.5771, 14.5821, 10.9659, 5.4828, -0.6240, -4.8980, -6.7653, -7.8702, -9.5130, -11.2193, -12.1947, -12.7875, -13.5035, -13.5705, -12.0467, -10.1946, -11.1638, -16.5918, -24.4423, -31.0464, -34.6433, -35.9300, -35.6345, -33.1863, -28.4086, -23.3184, -20.7883, -21.3552, -22.2727, -20.1604, -14.2286, -6.6712, -0.2897, 3.9291, 7.2435, 11.7306, 18.3403, 25.8358, 31.5584, 33.6963, 33.0131, 31.9688, 31.9311, 31.8425, 30.2950, 28.4476, 29.5386, 34.9635, 41.9133, 46.2495, 47.5332, 49.7347, 56.2082, 64.9135, 69.9009, 67.4579, 59.7465, 52.0974, 47.6860, 45.5773, 43.2243, 39.5108, 35.1334, 31.1330, 27.9661, 25.6310, 23.7599, 21.4133, 17.4881, 11.9040, 6.2104, 2.4992, 1.4770, 1.7193, 0.9566, -1.8908, -5.9543, -9.6137, -12.3488, -15.3076, -19.9959, -26.4781, -33.0153, -37.5012, -39.1330, -38.6878, -37.5101, -36.5794, -36.3929, -37.1123, -38.4073, -39.4113, -39.3437, -38.2809, -37.0347, -36.2082, -35.7034, -35.3572, -35.6986, -37.5373, -40.7496, -43.9003, -45.3575, -44.6201, -42.3697, -39.4889, -36.5009, -33.9639, -32.7754, -33.4650, -35.2459, -36.3714, -35.8743, -34.7958, -35.2138, -37.7878, -40.3913, -39.4429, -32.7844, -21.4514, -8.9438, 0.9904, 6.1089, 6.1292, 2.4378, -2.3013, -4.8182, -2.5870, 4.5964, 14.0697, 21.6238, 24.0321, 20.9710, 14.8400, 8.6489, 3.7782, -0.4566, -5.1463, -9.7612, -11.9396, -9.5785, -3.3168, 3.2273, 6.1762, 4.7136, 1.5217, -0.5612, -2.2666, -7.3886, -17.7591, -29.2673, -34.5008, -30.2512, -21.5223, -16.7207, -18.9577, -23.4166, -23.6684, -19.1384, -15.1181, -16.3891, -22.7523, -30.7493, -37.4560, -41.1877, -40.2336, -33.8174, -24.7970, -19.1177, -20.5891, -26.6461, -30.8713, -30.0455, -27.5928, -29.2882, -36.1963, -43.0645, -44.1063, -39.5927, -35.5159, -36.2770, -38.6955, -34.9969, -22.4044, -8.3644, -4.1846, -13.1775, -26.8115, -33.8706, -32.3924, -30.6639, -36.4345, -47.2267, -52.8368, -46.7304, -33.4708, -24.6068, -27.9011, -40.5155, -51.6275, -51.1121, -36.7782, -15.2168, 3.1098, 10.1852, 4.2453, -9.9225, -23.2828, -26.7960, -16.7758, 1.6186, 16.7396, 18.7063, 7.3183, -7.7101, -14.8110, -9.7809, 1.6305, 9.1867, 6.6682, -3.2355, -11.4136, -10.6400, -2.7631, 2.8362, -0.4647, -7.8795, -7.1957, 6.9212, 24.0017, 26.0517, 5.9168, -23.7385, -41.0058, -34.4850, -12.9289, 3.4885, 0.9503, -17.4524, -36.7747, -45.0977, -43.7659, -42.8587, -47.0129, -48.5099, -37.5461, -17.9555, -6.3646, -11.5005, -18.6721, -3.7817, 34.5674, 64.6413, 52.5834, 4.1036, -33.3991, -17.5697, 41.4199, 86.4835, 69.9434, 2.0114, -59.7048, -67.4966, -26.4956, 18.9470, 33.9745, 20.4952, -1.2782, -23.8018, -52.0198, -76.8150, -68.0323, -8.1607, 72.7581, 112.6802, 75.6408, -12.0744, -89.3151, -117.6612, -104.7067, -77.3307, -48.9987, -19.9650, 3.0824, 4.3859, -21.5597, -57.6276, -77.2484, -68.1581, -41.6151, -21.1683, -23.1948, -43.3733, -58.8777, -48.6828, -17.4014, 4.7655, -5.7213, -34.6026, -42.3076, -12.9083, 19.6414, 6.7807, -55.8015, -118.0950, -128.8001, -87.7729, -37.4437, -10.8411, -3.3353, 3.7364, 8.1759, -10.2889, -55.5504, -96.5260, -90.8639, -26.9755, 58.8598, 107.3746, 82.0146, -1.9464, -88.6477, -125.3889, -103.4100, -60.4166, -42.8351, -63.7707, -96.1267, -104.2293, -79.8251, -47.8772, -40.2633, -64.9444, -98.5198, -107.1824, -77.3312, -29.5853, -3.0226, -21.5528, -71.7234, -113.1799, -113.5511, -76.2119, -36.1186, -29.2159, -64.1456, -119.3274, -162.7172, -173.1820, -147.6990, -96.7829, -38.7146, 5.3801, 21.1957, 12.3445, -2.5991, -11.6661, -25.7949, -62.4172, -109.4214, -119.5496, -57.8159, 46.6795, 113.7883, 87.4756, -9.7322, -105.3489, -150.0506, -147.4447, -118.9106, -68.6380, -0.4732, 49.8546, 33.0040, -54.2900, -147.8416, -173.0396, -119.8138, -48.7928, -24.0836, -51.5651, -86.5522, -91.1058, -72.7276, -67.3982, -92.6137, -123.4519, -117.9475, -66.4827, -10.7384, -3.8774, -48.8339, -86.8775, -60.2503, 16.9410, 66.6040, 32.6408, -53.0810, -104.9302, -79.1996, -18.6014, 3.6999, -30.7485, -71.6930, -66.4665, -22.2520, 6.3320, -17.4146, -71.8869, -102.4869, -78.0262, -16.8595, 37.5409, 58.0019, 49.1563, 29.9789, 11.3301, -6.1926, -19.3006, -17.5459, 3.9195, 31.8886, 42.6582, 24.1772, -12.0430, -42.2603, -50.8979, -38.4411, -14.4927, 9.4583, 20.5858, 7.7709, -27.2928, -64.4212, -79.7452, -69.8820, -55.2788, -53.8373, -55.8279, -34.7891, 15.0934, 61.7649, 63.6061, 13.0707, -53.7965, -90.8339, -79.8123, -37.0597, 7.4552, 31.7053, 30.1192, 13.2970, 1.6024, 10.2032, 33.8757, 49.8568, 41.0535, 14.7475, -6.3973, -10.1520, -5.0786, -4.8358, -10.5007, -14.1414, -14.0559, -16.5655, -22.6047, -22.3281, -8.7720, 9.0957, 15.0725, 5.6965, -6.4932, -10.5103, -9.0506, -8.0158, -1.6920, 20.4046, 52.6474, 71.2079, 59.8582, 31.6149, 15.8381, 25.5942, 46.4460, 58.6292, 60.9421, 66.1536, 78.7175, 87.3932, 81.2183, 65.2597, 55.8765, 63.5033, 83.6583, 102.6978, 108.8094, 98.9333, 79.5469, 61.7355, 52.7166, 50.9209, 50.6153, 50.6067, 54.9512, 64.5547, 72.6305, 71.7191, 62.0721, 48.9671, 34.4707, 18.2849, 6.1948, 8.9152, 27.2880, 44.4685, 41.4769, 20.0274, 2.3746, 6.1090, 25.0257, 38.9213, 37.3145, 26.9482, 18.0800, 11.6734, 3.2528, -5.8838, -7.1206, 5.4528, 25.5412, 36.2430, 23.4701, -7.9509, -31.5367, -20.1499, 23.3666, 60.6462, 53.1858, 5.0228, -38.9314, -41.7122, -12.5081, 12.7492, 20.8530, 33.8497, 69.1672, 104.6393, 99.2191, 44.2588, -20.3915, -44.8886, -15.9543, 36.7996, 71.8023, 68.6003, 40.4760, 24.5749, 53.4220, 121.7893, 180.2999, 175.3698, 104.1724, 23.0936, -7.1464, 19.0476, 55.1212, 59.5197, 35.9743, 15.1022, 12.0066, 15.8652, 14.7804, 13.1557, 18.9569, 25.5119, 17.4138, -7.8287, -32.9693, -38.0035, -19.5692, 7.9439, 26.7838, 29.6653, 23.2559, 21.3644, 32.2589, 48.7502, 52.1240, 30.8832, -3.2005, -21.9441, -7.1845, 30.7069, 62.5225, 66.7678, 44.8800, 14.0574, -8.0199, -11.0520, 5.8516, 31.1000, 42.1800, 21.4008, -23.6836, -60.5359, -60.2218, -26.1817, 10.8941, 27.7725, 31.5260, 42.8304, 62.3740, 65.6605, 34.9758, -12.6349, -39.7214, -29.0481, 0.8129, 23.7154, 38.0893, 58.3581, 79.7151, 72.2475, 22.2220, -32.9702, -33.1330, 35.7824, 115.4537, 130.8231, 67.8278, -14.7254, -49.2242, -24.3889, 15.0193, 21.3695, -11.9200, -53.2624, -68.3227, -51.5038, -28.7502, -30.5770, -58.4340, -79.1406, -59.4610, -6.0927, 38.0451, 40.7358, 13.1618, -12.5512, -27.3825, -48.3382, -78.4580, -87.1222, -46.6575, 21.9124, 60.5443, 35.4049, -26.2138, -66.4936, -56.3736, -17.6997, 6.7180, -4.8532, -38.8988, -64.3817, -61.0734, -35.6873, -14.1685, -16.6226, -38.5592, -58.3181, -61.9038, -55.5481, -52.8400, -56.4495, -58.4606, -54.4145, -48.2209, -42.7070, -33.5860, -17.6921, -2.6862, -2.4278, -22.6678, -53.9881, -78.1808, -80.1648, -57.3206, -22.9633, -0.0864, -4.5168, -29.7902, -51.6337, -51.2644, -32.7683, -14.9147, -9.2100, -12.0884, -17.6411, -28.8489, -50.0201, -72.9057, -79.8882, -63.6159, -37.6974, -22.9948, -25.7812, -34.6955, -38.7634, -41.6936, -53.8967, -74.2335, -87.0994, -78.9765, -52.0529, -18.1618, 14.4539, 40.0971, 45.1356, 13.9180, -45.5631, -91.2198, -79.1325, -10.6728, 62.8532, 89.5670, 64.7568, 23.4692, -4.9667, -19.7176, -27.6985, -20.7676, 11.4198, 53.8313, 69.4864, 36.3400, -26.5871, -74.3803, -74.0756, -27.6860, 33.2352, 69.9115, 60.5710, 14.4435, -35.8333, -61.2011, -57.7660, -40.0953, -17.9195, 10.0832, 38.3607, 44.8062, 12.8509, -40.2658, -71.8503, -56.4034, -12.8383, 16.7947, 10.9729, -15.5144, -36.0142, -42.0341, -44.2810, -49.1858, -43.4040, -7.4049, 55.2316, 107.4046, 106.1919, 46.5219, -26.3786, -56.1747, -28.5604, 20.6016, 51.0024, 60.7336, 73.5239, 91.2343, 80.5714, 18.4873, -61.0480, -86.4298, -22.1306, 85.4132, 145.0451, 104.9507, -1.7689, -91.1744, -102.6025, -42.0176, 36.5482, 82.8639, 84.8992, 61.8175, 34.1802, 7.0203, -20.1374, -36.6491, -25.3090, 13.6373, 51.3728, 54.2833, 20.7422, -13.6315, -13.8036, 15.0234, 33.4223, 14.2889, -24.3912, -39.6408, -13.1171, 27.6643, 40.4712, 14.9382, -18.5452, -25.6619, -6.6087, 7.6448, -5.1093, -30.6699, -35.5258, -8.2124, 24.9378, 29.4080, 2.4873, -24.5527, -24.5882, -6.4541, -2.2647, -27.2261, -61.3363, -73.1062, -54.4331, -26.6120, -15.5949, -27.6637, -48.2269, -57.3198, -44.0561, -12.4410, 20.5902, 35.0894, 23.1526, -3.3953, -23.6998, -28.3192, -26.1924, -28.8393, -31.3323, -16.0228, 22.4666, 61.4777, 67.4083, 31.2281, -17.4374, -37.5068, -18.0868, 12.1388, 19.8428, 5.7290, -0.0853, 19.6918, 44.0904, 39.7610, 5.6467, -22.7353, -16.4827, 11.9656, 26.8778, 13.6530, -9.4092, -22.4141, -29.7717, -45.4187, -64.0553, -64.2186, -39.8491, -13.9436, -13.1228, -36.7598, -60.2904, -64.9912, -54.8457, -44.3280, -39.3609, -35.3891, -28.8017, -23.7437, -27.9290, -44.1365, -64.3603, -71.3942, -51.9265, -14.6921, 7.9482, -11.8521, -62.9909, -97.7980, -78.7383, -20.6287, 22.7890, 15.3701, -26.6005, -58.7931, -57.9776, -37.1298, -17.5048, -1.8585, 18.2043, 39.4572, 45.8683, 30.9745, 11.1984, 10.1138, 32.3900, 58.7009, 66.7746, 54.5123, 41.1323, 46.3918, 69.2204, 87.3554, 79.7474, 49.9649, 24.5201, 24.7931, 42.7181, 50.1937, 33.7859, 13.3943, 20.1249, 58.2004, 96.6831, 101.4514, 70.2497, 32.1768, 15.9432, 23.8705, 36.0657, 34.1640, 18.2634, 3.7917, 5.2672, 22.4353, 40.0298, 41.2233, 22.9171, -1.8135, -16.5646, -15.9281, -6.8180, 1.3757, 4.2298, 1.2007, -6.9314, -16.3676, -19.9943, -14.3686, -6.6507, -9.4896, -26.8509, -48.1022, -58.7742, -54.8766, -43.5349, -31.6647, -19.3527, -5.4515, 5.3612, 5.7307, -5.2459, -19.7329, -28.5887, -28.7920, -23.1009, -15.2774, -7.8171, -2.6717, -0.9744, -1.3561, -0.3324, 4.5383, 12.5726, 21.4050, 29.3542, 35.4910, 38.5476, 37.7001, 34.5170, 32.3038, 32.3337, 31.7152, 26.6293, 17.6355, 10.0420, 7.8312, 8.1281, 3.7111, -8.1753, -20.7935, -23.5272, -12.9180, 3.0063, 12.1775, 9.6000, 0.3526, -7.3751, -10.1587, -10.0744, -10.0370, -10.6211, -10.8069, -9.6259, -6.6435, -2.4737, 0.5652, 0.0482, -3.0334, -3.9089, 0.7462, 8.1911, 12.9474, 14.5750, 18.5909, 28.0101, 37.4478, 39.7074, 36.7155, 38.3889, 48.7535, 57.7123, 51.3681, 29.3997, 7.6582, 2.2963, 14.2421, 30.8599, 40.0997, 40.0386, 36.0907, 33.1226, 31.8309, 30.4279, 27.4194, 22.9402, 18.5415, 15.6141, 13.6046, 10.5632, 6.2950, 3.9223, 6.4363, 11.9030, 14.5086, 11.6697, 8.3737, 12.0985, 23.6649, 35.4222, 39.5677, 36.9237, 35.9485, 43.3605, 56.9625, 67.9417, 69.4573, 62.6629, 54.9221, 52.6494, 55.6411, 57.8973, 53.8194, 43.6083, 32.6311, 25.9345, 24.0964, 24.3546, 24.3630, 23.5767, 21.7622, 18.4166, 14.5408, 13.2311, 16.3546, 20.8560, 20.4909, 12.3904, 0.9089, -6.3516, -6.1215, -2.1865, -0.8137, -4.6224, -11.1249, -15.9336, -16.2969, -12.2263, -5.8632, -0.3423, 1.7350, -0.0458, -3.4895, -5.4895, -4.6503, -2.5974, -2.9674, -8.9514, -21.0024, -36.0499, -48.9401, -55.5995, -55.8061, -53.1391, -51.8856, -53.5847, -56.2408, -56.5312, -52.6631, -45.6177, -38.4590, -34.5473, -35.4118, -39.2300, -41.4422, -38.0182, -28.9935, -18.5597, -11.0743, -7.0603, -3.5713, 1.8311, 8.5362, 14.6630, 20.2714, 27.4731, 37.3866, 48.1380, 56.4914, 60.8803, 62.2442, 62.3335, 62.0069, 61.1202, 59.2922, 56.4660, 53.1420, 50.2298, 48.1579, 45.8667, 41.3749, 34.1280, 26.5132, 22.1356, 22.3959, 25.4253, 28.6221, 31.4761, 35.3033, 40.8574, 47.2897, 53.4545, 59.2175, 65.0100, 70.8202, 76.3597, 81.8322, 87.6312, 93.2127, 97.1019, 98.5682, 98.9080, 100.3866, 103.7182, 106.8466, 106.4905, 101.1520, 92.9238, 86.3819, 85.1492, 88.7983, 92.9778, 93.1014, 88.3648, 82.0959, 78.0921, 76.9637, 76.3499, 74.1948, 70.8602, 67.8938, 65.7670, 63.7982, 61.9449, 61.4381, 63.0905, 65.6353, 66.5389, 64.3826, 59.9232, 54.8831, 50.3580, 46.6453, 43.9664, 42.5699, 42.1233, 41.5793, 39.9580, 37.0768, 33.3407, 29.0092, 23.9641, 18.2563, 12.7473, 9.0316, 8.4348, 10.7100, 13.6051, 14.1229, 10.8211, 5.1312, 0.2237, -1.6823, -1.0051, 0.0197, -0.3994, -2.1325, -3.4296, -2.2378, 2.5913, 10.6076, 19.6265, 26.4711, 28.7563, 26.7841, 23.7768, 23.6058, 27.5661, 33.1688, 36.3463, 35.0655, 30.7958, 26.5805, 24.0822, 22.6947, 21.0262, 18.5299, 15.6903, 13.3387, 12.1149, 11.9620, 11.4840, 8.2451, 0.8346, -8.9925, -16.9803, -19.9155, -18.4320, -16.0435, -15.7728, -18.1879, -22.1656, -26.2803, -29.0639, -29.0925, -26.1389, -22.1113, -19.7296, -19.8223, -20.5658, -19.8258, -17.6095, -15.4079, -13.4814, -10.0061, -3.4392, 4.9563, 11.7140, 14.1111, 11.9784, 6.8657, 0.4418, -5.8409, -10.3820, -12.0116, -11.2968, -10.6819, -12.7467, -18.1844, -25.3386, -31.2893, -33.2723, -29.7888, -21.4478, -11.1324, -2.7308, 1.2180, 1.1759, 0.0572, 0.9058, 4.6726, 9.6881, 13.1318, 13.4888, 12.0985, 12.1588, 15.6358, 21.0518, 24.8708, 25.3881, 24.7266, 26.2922, 30.4216, 33.4332, 31.4331, 24.5815, 17.0001, 12.6192, 11.6705, 11.2347, 8.4850, 2.8906, -4.2763, -11.6313, -18.5259, -24.7806, -30.2980, -35.0885, -39.3517, -43.2992, -46.9414, -50.0981, -52.4967, -53.7912, -53.6671, -52.1338, -49.6477, -46.7201, -43.3860, -39.2771, -34.2854, -28.9413, -23.9792, -19.6747, -15.7366, -11.6378, -6.8393, -1.0238, 5.2775, 10.3270, 12.3193, 11.5215, 10.7843, 12.9151, 17.3803, 20.4733, 19.3723, 15.5423, 13.2833, 14.9440, 18.3783, 19.4985, 16.7090, 12.0843, 8.3901, 5.9390, 2.9702, -1.3889, -5.6177, -7.7284, -8.1266, -9.5693, -13.8917, -19.3434, -22.0207, -20.1484, -16.4030, -15.4896, -19.5395, -26.2768, -31.6259, -33.5668, -33.1854, -32.4707, -31.8662, -30.1958, -26.5226, -21.5350, -16.9940, -14.0323, -12.1570, -9.8523, -5.9928, -0.6348, 5.3426, 10.9875, 15.4759, 17.8372, 17.0892, 13.1493, 7.6724, 3.6163, 3.4314, 7.2208, 12.5016, 15.9395, 15.7089, 12.5382, 8.6782, 5.9395, 4.6302, 4.0045, 3.2241, 1.7488, -0.7318, -4.1848, -8.0993, -11.9788, -16.1068, -21.2903, -27.4145, -32.6682, -35.0105, -34.5450, -33.7921, -35.1566, -38.4779, -41.5626, -42.8605, -42.9080, -43.0908, -43.8050, -44.4288, -44.6373, -44.8542, -45.3140, -45.4408, -44.6544, -43.3764, -42.6060, -42.5884, -42.4988, -41.5631, -39.9234, -38.0332, -35.7044, -32.5043, -29.0536, -27.1229, -27.8988, -30.4302, -32.2956, -31.8383, -29.5123, -27.1754, -26.4931, -28.0302, -31.2728, -35.0018, -37.8052, -38.7577, -37.8459, -35.7057, -33.0672, -30.6394, -29.2428, -29.3243, -30.0151, -29.2364, -25.5347, -19.9079, -15.2771, -13.8288, -15.0872, -16.7630, -17.1688, -16.4565, -15.7463, -15.6522, -15.8238, -15.5265, -14.3911, -12.8487, -12.0887, -13.2617, -16.1800, -18.8799, -19.1841, -17.1421, -15.5316, -17.2921, -22.3703, -27.4788, -29.0875, -26.3990, -21.4342, -16.8446, -14.1105, -13.3897, -14.1231, -15.4442, -16.5665, -17.3021, -18.0008, -18.6249, -18.1759, -15.6585, -11.6249, -8.1364, -6.8681, -7.6148, -8.9353, -9.7944, -9.9937, -9.3037, -7.1103, -3.2757, 1.2036, 4.8799, 6.9478, 7.3493, 6.2192, 3.9256, 1.6003, 0.8405, 2.4190, 5.5459, 8.6346, 10.5048, 10.6915, 9.2643, 7.1860, 6.5952, 9.4925, 15.3959, 20.8515, 22.2382, 19.2310, 15.0034, 12.8668, 13.1322, 13.4442, 11.6293, 7.5864, 2.5652, -2.5978, -8.0320, -13.9197, -19.6233, -24.0216, -26.6894, -28.6291, -31.6453, -36.6573, -42.3559, -45.6542, -43.9761, -37.5350, -29.4009, -23.1152, -20.0735, -19.0915, -18.1801, -16.3456, -13.7122, -10.5743, -7.0726, -3.8423, -2.2856, -3.5076, -6.9076, -10.2443, -11.4576, -10.4542, -8.9665, -8.6490, -9.4549, -9.9858, -9.3309, -8.2770, -8.4525, -10.2811, -12.0614, -11.3878, -7.6054, -2.8084, -0.3471, -2.2629, -7.9147, -14.6894, -19.8419, -21.9582, -21.2949, -19.2008, -17.1983, -16.1789, -15.9390, -15.1982, -12.2608, -6.2375, 1.9159, 9.6966, 14.9712, 17.9735, 21.1800, 26.7746, 34.1895, 40.4437, 43.0228, 42.2033, 40.4329, 39.7192, 39.9513, 39.6164, 37.6304, 34.1993, 30.1831, 26.0956, 21.8718, 17.3825, 12.9098, 9.0692, 6.3555, 4.7720, 3.7680, 2.4332, -0.1145, -4.1850, -9.0921, -13.3225, -15.4070, -14.8419, -12.2184, -8.5053, -4.3016, 0.2524, 5.0767, 9.8568, 14.0760, 17.2197, 19.0630, 19.9730, 20.8764, 22.6128, 25.1411, 27.4751, 28.4805, 27.6799, 25.4179, 22.6946, 21.0127, 21.7652, 24.9201, 28.3341, 29.3137, 27.4761, 25.5302, 26.3320, 29.2089, 30.1005, 25.8179, 17.6551, 10.1947, 6.7948, 6.8691, 7.4784, 6.6900, 4.8606, 3.2293, 2.2604, 1.6140, 1.1092, 1.0740, 1.7196, 2.5312, 2.5572, 1.3950, -0.1132, -0.3466, 1.5620, 4.4582, 5.7384, 3.7204, -0.4782, -3.8952, -4.8629, -4.8502, -6.6614, -10.9216, -15.0146, -15.7906, -12.8508, -8.6707, -5.7940, -4.6344, -4.1200, -3.6977, -3.9619, -5.4650, -7.6433, -9.3451, -10.2363, -11.3158, -13.8729, -17.9455, -21.8837, -23.5854, -22.3150, -19.3013, -16.4393, -14.4546, -12.5672, -9.9834, -7.3391, -6.2549, -7.4686, -9.7365, -10.7849, -9.3221, -6.1629, -3.6002, -3.7451, -7.1589, -12.6559, -18.2810, -22.6652, -25.7096, -28.1867, -30.8333, -33.7980, -36.6607, -38.6630, -38.9712, -37.1317, -33.5657, -29.4851, -26.0164, -23.3140, -20.6306, -17.2581, -13.2646, -9.2787, -5.7029, -2.2754, 1.6235, 6.2921, 11.2776, 15.5471, 18.1931, 19.2266, 19.8531, 21.8797, 26.5053, 33.2649, 39.9145, 43.5061, 42.0878, 35.9004, 27.1612, 18.5920, 11.8878, 7.2652, 4.0919, 1.6806, -0.3981, -2.2916, -3.7654, -4.2283, -3.2478, -1.2369, 0.5521, 0.9843, 0.0511, -1.1730, -1.6196, -1.2125, -0.7313, -0.7428, -0.7432, 0.5812, 4.1578, 9.4834, 14.7387, 18.0672, 18.9592, 18.5273, 18.4521, 19.6070, 21.6108, 23.5327, 24.8284, 25.5542, 25.8922, 25.7212, 24.6432, 22.2185, 18.1503, 12.5634, 6.3787, 1.1999, -1.6974, -2.5062, -2.7988, -4.0791, -6.5397, -9.3740, -12.0010, -14.5294, -17.1325, -19.6085, -21.8232, -24.0000, -25.9131, -26.0650, -22.6537, -15.8861, -8.8641, -5.4011, -6.7046, -10.3896, -12.8063, -12.1979, -9.6595, -7.3838, -6.2288, -4.9154, -1.4526, 4.6212, 11.3406, 15.7504, 16.3716, 14.2074, 11.2859, 8.5260, 5.4697, 2.0314, -0.2585, 0.5035, 4.1173, 7.6820, 8.0004, 4.5435, -0.2885, -3.9178, -6.3250, -9.8353, -16.3352, -25.0405, -33.1119, -38.2097, -40.2084, -40.6357, -40.9925, -41.8815, -43.3489, -45.3534, -47.6202, -49.4235, -50.0760, -49.7930, -49.8144, -51.3885, -54.5828, -58.1824, -60.7873, -61.9862, -62.4896, -63.2355, -64.4321, -65.3964, -65.2078, -63.5170, -60.8281, -58.0458, -55.6846, -53.4846, -50.8811, -47.8608, -45.1631, -43.4628, -42.5207, -41.4763, -39.9658, -38.5493, -37.8911, -37.9885, -38.6562, -40.3957, -43.9319, -48.6783, -52.3498, -52.6701, -49.2131, -43.1580, -35.5982, -26.8097, -17.0843 };
	char *name = "db4";
	obj = wave_init(name);// Initialize the wavelet
	N = 512;  //Length of Signal
	inp = (float*)malloc(sizeof(float)* N); //Input signal
	out = (float*)malloc(sizeof(float)* N); //Output signal

	//printf("Input vector >>>>>>>>>>>>>\n");
	for (i = 0; i < N; ++i) {
		inp[i] = temp[i];
	}

	J = 8; //Decomposition Levels
	
	wt = wt_init(obj, N, J);// Initialize the wavelet transform object

	for (i = 0; i < 60*480; ++i) {
		dwt(wt, inp);// Perform DWT
	}
	// DWT output coefficients can be accessed using wt->output vector
	for (i = 0; i < wt->outlength; ++i) {
			wt->output[i] = 50.0f*(int)(wt->output[i] / 50.0f);
	}

	//for (i = 0; i < 10 * 60 * 480; ++i) {
		idwt_sym_direct(wt, out);// Perform IDWT
	//}

	// Test reconstruction output
	printf("Input	Output	Wavelet\n");
	for (i = 0; i < N; ++i) {
		printf("%g	", inp[i]);
		printf("%g	", out[i]);
		printf("%g\n", wt->output[i]);
	}

	wave_free(obj);
	wt_free(wt);
	free(inp);
	free(out);
	return 0;
}