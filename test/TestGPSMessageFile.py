from gpsMessageFile import GpsMessageFile
from gpsMessageFile import GPSFactory
from gpsMessageFile import RinexFileReader
from nose.tools import assert_equals

fileContent = ["G01 2017 02 13 22 00 00 4.918128252029E-05 9.094947017729E-13 0.000000000000E+00",
               "     6.000000000000E+00-7.018750000000E+01 4.509830709694E-09-2.275517339128E+00",
               "    -3.596767783165E-06 6.229116464965E-03 4.574656486511E-06 5.153692497253E+03",
               "     1.656000000000E+05 9.313225746155E-09-6.361854098473E-01-5.960464477539E-08",
               "     9.666622276019E-01 2.948750000000E+02 5.418360356115E-01-8.197841473003E-09",
               "    -2.078658013021E-10 1.000000000000E+00 1.936000000000E+03 0.000000000000E+00",
               "     2.000000000000E+00 0.000000000000E+00 5.122274160385E-09 6.000000000000E+00",
               "     1.584180000000E+05 4.000000000000E+00"]
gps = GPSFactory.createGPSFromRinexFile(fileContent,3)
fileContent2 = [' 1 02  1 10  0  0  0.0 2.108854241669E-04 1.591615728103E-12 0.000000000000E+00',
     '    1.740000000000E+02-6.625000000000E+01 3.994452099249E-09 2.453913345123E+00',
     '   -3.479421138763E-06 5.274268332869E-03 1.066364347935E-05 5.153708732605E+03',
     '    3.456000000000E+05 1.378357410431E-07-2.987566298300E-01-2.235174179077E-08',
     '    9.676350812908E-01 1.789375000000E+02-1.727933015620E+00-7.573886911362E-09',
     '    3.507288949806E-10 0.000000000000E+00 1.148000000000E+03 0.000000000000E+00',
     '    2.800000000000E+00 0.000000000000E+00-3.259629011154E-09 4.300000000000E+02',
     '    3.455400000000E+05']

gps2 = GPSFactory.createGPSFromRinexFile(fileContent2,2)


fileContentReceptor = ["  5 11  3 16  0  0.0000000  0  8G14G19G 1G 3G22G20G11G25",
"  21286668.455   111862187.04947  21286668.514    87165341.44247",
"  21267496.922   111761439.40447  21267495.926    87086832.71248",
"  20081605.043   105529566.39347  20081606.065    82230836.41048",
"  22593501.765   118729615.82047  22593505.077    92516581.20747",
"  23663167.973   124350739.80947  23663170.114    96896675.78146",
"  24096330.954   126627077.26446  24096336.220    98670478.13946",
"  22054493.203   115897167.96247  22054494.097    90309499.84847",
"  22879205.847   120231083.78046  22879210.454    93686594.71447"]

receptor = GPSFactory.createReceptorFromRinexFile(fileContentReceptor, (1,2,3), {' 2':14000},2)

def test_13():
    r = RinexFileReader()
    receptores = r.readReceptorFromObservation('test/UNISINOS/Base3070.05o')
    satelites = r.readGPSFromEphemeris('test/UNISINOS/Base3070.05n')
    #print('Receptor',receptores[0].sat_number,receptores[0].epochYear,receptores[0].epochMonth,receptores[0].epochDay,receptores[0].epochHour,receptores[0].epochMinute,receptores[0].epochSecond)
    #for i in range(len(satelites)):
    #    print(satelites[i].sat_number,satelites[i].epochYear,satelites[i].epochMonth,satelites[i].epochDay,satelites[i].epochHour,satelites[i].epochMinute,satelites[i].epochSecond)
    dic = r.transformGPSListDictionary(satelites)
    receptores[0].loadSateliteData(dic)
    #print(receptores[0].getCoordinates())
    assert_equals(len(receptores[0].gps),8)


def test_12():
    assert_equals(receptor.epochYear, 5)
    assert_equals(receptor.epochMonth, 11)
    assert_equals(receptor.epochDay, 3)
    assert_equals(receptor.epochHour, 16)
    assert_equals(receptor.epochMinute, 0)
    assert_equals(receptor.epochSecond, 0)

def test_11():
    assert_equals(receptor.aprox_pos, (1,2,3))

def test_10():
    gps2.coord_mult = -3
    gps2.dtr = 900
    assert_equals(gps2.coordinate_WGS84(), (20573.622159228864,   3151.6233033090466,  16691.438710244893))

def test_09():
    gps2.coord_mult = -3
    gps2.dtr = 0
    assert_equals(gps2.coordinate_WGS84(), (22137.659080861868, 2318.121737218721, 14690.038045001766))


#eighth line
def test_08():
    assert_equals(gps.transTime, 1.584180000000E+05)
    assert_equals(gps.fitInterval, 4.000000000000E+00)


#seventh line
def test_07():
    assert_equals(gps.svAccuracy, 2.000000000000E+00)
    assert_equals(gps.svHealth, 0.000000000000E+00)
    assert_equals(gps.tgd, 5.122274160385E-09)
    assert_equals(gps.iodc, 6.000000000000E+00)

#sixth line
def test_06():
    assert_equals(gps.idot, -2.078658013021E-10)
    assert_equals(gps.codesL2Channel, 1.000000000000E+00)
    assert_equals(gps.gpsWeek, 1.936000000000E+03)
    assert_equals(gps.l2PDataFlag, 0.000000000000E+00)

#fifth line
def test_05():
    assert_equals(gps.i0, 9.666622276019E-01)
    assert_equals(gps.crc, 2.948750000000E+02)
    assert_equals(gps.omega, 5.418360356115E-01)
    assert_equals(gps.omegaDot, -8.197841473003E-09)

#fourth line
def test_04():
    assert_equals(gps.toeTime, 1.656000000000E+05)
    assert_equals(gps.cic, 9.313225746155E-09)
    assert_equals(gps.omega0, -6.361854098473E-01)
    assert_equals(gps.cis, -5.960464477539E-08)

#third line
def test_03():
    assert_equals(gps.cuc, -3.596767783165E-06)
    assert_equals(gps.eccentricity, 6.229116464965E-03)
    assert_equals(gps.cus, 4.574656486511E-06)
    assert_equals(gps.sqrtA, 5.153692497253E+03)

#second line
def test_02():
    assert_equals(gps.iode, 6.000000000000E+00)
    assert_equals(gps.crs, -7.018750000000E+01)
    assert_equals(gps.deltaN, 4.509830709694E-09)
    assert_equals(gps.m0, -2.275517339128E+00)

#first line
def test_01():
    assert_equals(gps.sat_number, "G01")
    assert_equals(gps.epochYear, 2017)
    assert_equals(gps.epochMonth, 2)
    assert_equals(gps.epochDay, 13)
    assert_equals(gps.epochHour, 22)
    assert_equals(gps.epochMinute, 0)
    assert_equals(gps.epochSecond, 0)
    assert_equals(gps.sv_clock_bias, 4.918128252029E-05)
    assert_equals(gps.sv_clock_drift, 9.094947017729E-13)
    assert_equals(gps.sv_clock_drift_rate, 0.000000000000E+00)
