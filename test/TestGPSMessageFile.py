from gpsMessageFile import GpsMessageFile
from gpsMessageFile import GPSFactory
from nose.tools import assert_equals

fileContent = ["G01 2017 02 13 22 00 00 4.918128252029E-05 9.094947017729E-13 0.000000000000E+00",
               "     6.000000000000E+00-7.018750000000E+01 4.509830709694E-09-2.275517339128E+00",
               "    -3.596767783165E-06 6.229116464965E-03 4.574656486511E-06 5.153692497253E+03",
               "     1.656000000000E+05 9.313225746155E-09-6.361854098473E-01-5.960464477539E-08",
               "     9.666622276019E-01 2.948750000000E+02 5.418360356115E-01-8.197841473003E-09",
               "    -2.078658013021E-10 1.000000000000E+00 1.936000000000E+03 0.000000000000E+00",
               "     2.000000000000E+00 0.000000000000E+00 5.122274160385E-09 6.000000000000E+00",
               "     1.584180000000E+05 4.000000000000E+00"]
gps = GPSFactory.createGPSFromRinexFile(fileContent,'3')
fileContent2 = [' 1 02  1 10  0  0  0.0 2.108854241669E-04 1.591615728103E-12 0.000000000000E+00',
     '    1.740000000000E+02-6.625000000000E+01 3.994452099249E-09 2.453913345123E+00',
     '   -3.479421138763E-06 5.274268332869E-03 1.066364347935E-05 5.153708732605E+03',
     '    3.456000000000E+05 1.378357410431E-07-2.987566298300E-01-2.235174179077E-08',
     '    9.676350812908E-01 1.789375000000E+02-1.727933015620E+00-7.573886911362E-09',
     '    3.507288949806E-10 0.000000000000E+00 1.148000000000E+03 0.000000000000E+00',
     '    2.800000000000E+00 0.000000000000E+00-3.259629011154E-09 4.300000000000E+02',
     '    3.455400000000E+05']
gps2 = GPSFactory.createGPSFromRinexFile(fileContent2,'2')

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
    assert_equals(gps.epochYear, "2017")
    assert_equals(gps.epochMonth, "02")
    assert_equals(gps.epochDay, "13")
    assert_equals(gps.epochHour, "22")
    assert_equals(gps.epochMinute, "00")
    assert_equals(gps.epochSecond, "00")
    assert_equals(gps.sv_clock_bias, 4.918128252029E-05)
    assert_equals(gps.sv_clock_drift, 9.094947017729E-13)
    assert_equals(gps.sv_clock_drift_rate, 0.000000000000E+00)
