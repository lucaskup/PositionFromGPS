import math
class RinexFileReader:
    def readGPSFromEphemeris(self,fileType,fileName):
        l = []
        with open(fileName,'r') as file:
            content = file.readlines()
            file.seek(0)
            endOfHeaderLineNumber = 0
            for i in range(len(content)):
                if 'END OF HEADER' in content[i]:
                    endOfHeaderLineNumber = i
                    break
            iterationNumber = (len(content) - endOfHeaderLineNumber-1) // 8
            for i in range(iterationNumber):
                posToCut = endOfHeaderLineNumber+1+(8*i)
                gpsData = content[posToCut:posToCut+9]
                for i in range(len(gpsData)):
                    gpsData[i] = gpsData[i].replace('D', 'E')
                l.append(GPSFactory.createGPSFromRinexFile(gpsData,fileType))
        return l

class GpsMessageFile:
    def deltaTsv(self):
        #PRECISA DE REVISAO
        #%tr   = (dia*24+hora)*3600+(minu*60)+seg;
        #tr=(0*24+0)*3600+(0*60)+30
        tr = self.toeTime
        #PD = 0
        #dtr = 0
        c = 299792458
        tgps = tr - self.PD/c;

        dts  = self.sv_clock_bias + self.sv_clock_drift*(tgps - self.toeTime) + self.sv_clock_drift_rate*(tgps-self.toeTime)**2;

        #%tempo de propagação aproximado do sinal
        tal = self.PD/c - self.dtr + dts

        tgps = tr - tal + dts;

        #%time diff from toe

        dt = tgps - self.toeTime;
        return dt
    def GM(self):
        return 3.986005E+14
    def omegaEarth(self):
        return 7.292115E-5
    def n0(self):
        return (self.GM()/(self.sqrtA**6))**0.5
    def n(self):
        return self.n0() + self.deltaN
    def M(self):
        return self.m0 + self.n()*self.deltaTsv()
    #Keplers equation through iteration
    def E(self):
        e0 = self.M()
        e1 = self.M() + self.eccentricity*math.sin(e0)
        lim = 0
        while e0 != e1 or lim >10:
            e0 = e1
            e1 = self.M() + self.eccentricity*math.sin(e0)
            lim += 1
        return e1
    #True Anomaly
    def V(self):
        #print('eccentricity',str(self.eccentricity))
        #print('E',str(self.E()))
        #print(str(math.acos((math.cos(self.E())-self.eccentricity)/(1-self.eccentricity*math.cos(self.E())))))
        #print(str(math.asin((((1-self.eccentricity**2)**0.5)*math.sin(self.E()))/(1-self.eccentricity*math.cos(self.E())))))
        #print(str(math.atan((((1-self.eccentricity**2)**0.5)*math.sin(self.E()))/(math.cos(self.E())-self.eccentricity))))
        return math.acos((math.cos(self.E())-self.eccentricity)/(1-self.eccentricity*math.cos(self.E())))
        #return math.acos((math.cos(self.E())-self.eccentricity)/(1-self.eccentricity*math.cos(self.E())))
    def VplusOmega(self):
        return self.V() + self.omega

    def tetaU(self):
        return self.cuc*math.cos(2*self.VplusOmega()) + self.cus*math.sin(2*self.VplusOmega())
    def U(self):
        return self.VplusOmega() + self.tetaU()

    def tetaR(self):
        return self.crc*math.cos(2*self.VplusOmega()) + self.crs*math.sin(2*self.VplusOmega())
    def R(self):
        return (self.sqrtA**2)*(1-self.eccentricity*math.cos(self.E())) + self.tetaR()

    def x_OrbitalPlane(self):
        return self.R()*math.cos(self.U())*10**self.coord_mult

    def y_OrbitalPlane(self):
        return self.R()*math.sin(self.U())*10**self.coord_mult

    def coordinate_OrbitalPlane(self):
        return (self.x_OrbitalPlane(),self.y_OrbitalPlane())

    def tetaI(self):
        return self.cic*math.cos(2*self.VplusOmega())+ self.cis*math.sin(2*self.VplusOmega())
    def I(self):
        return self.i0 + self.idot*self.deltaTsv() + self.tetaI()
    def capitalOmega(self):
        return self.omega0 + (self.omegaDot - self.omegaEarth())*self.deltaTsv()-self.omegaEarth()*self.toeTime

    def x_WGS84(self):
        return self.x_OrbitalPlane()*math.cos(self.capitalOmega()) - self.y_OrbitalPlane()*math.sin(self.capitalOmega())*math.cos(self.I())
    def y_WGS84(self):
        return self.x_OrbitalPlane()*math.sin(self.capitalOmega()) + self.y_OrbitalPlane()*math.cos(self.capitalOmega())*math.cos(self.I())
    def z_WGS84(self):
        return self.y_OrbitalPlane()*math.sin(self.I())
    def coordinate_WGS84(self):
        return (self.x_WGS84(),self.y_WGS84(),self.z_WGS84())
    def __init__(self):
        """dtr e o tempo que passou desde o horario importado que consta no arquivo importado"""
        self.dtr = 0
        """PD e a informacao de pseudo distancia que deve ser passada para tornar os calculos de posicao mais precisos"""
        self.PD = 0
        """o atributo coord_mult serve para aplicar um multiplicador 10**coord_mult ao resultado, permitindo variar
        a unidade da resposta entre metros e quilometros, por exemplo"""
        self.coord_mult = 0


"""Classe responsavel por criar instancias de GpsMessageFile de acordo com o Layout Rinex 2 ou Rinex 3"""
class GPSFactory:
    def createGPSFromRinexFile(fileBlock, version):
        if version == '2':
            g = GPSFactory._gpsFromRinex2(fileBlock)
        else:
            g = GPSFactory._gpsFromRinex3(fileBlock)
        return g
    def _gpsFromRinex2(fileBlock):
        g = GpsMessageFile()
        #first line
        g.sat_number = fileBlock[0][0:3]
        g.epochYear = fileBlock[0][4:6]
        g.epochMonth = fileBlock[0][7:9]
        g.epochDay = fileBlock[0][10:12]
        g.epochHour = fileBlock[0][13:15]
        g.epochMinute = fileBlock[0][16:18]
        g.epochSecond = fileBlock[0][19:23]
        g.sv_clock_bias = float(fileBlock[0][22:41])
        g.sv_clock_drift = float(fileBlock[0][41:60])
        g.sv_clock_drift_rate = float(fileBlock[0][60:79])

        #second line
        g.iode = float(fileBlock[1][0:22])
        g.crs = float(fileBlock[1][22:41])
        g.deltaN = float(fileBlock[1][41:60])
        g.m0 = float(fileBlock[1][60:79])

        #third line
        g.cuc = float(fileBlock[2][0:22])
        g.eccentricity = float(fileBlock[2][22:41])
        g.cus = float(fileBlock[2][41:60])
        g.sqrtA = float(fileBlock[2][60:79])

        #fourth line
        g.toeTime = float(fileBlock[3][0:22])
        g.cic = float(fileBlock[3][22:41])
        g.omega0 = float(fileBlock[3][41:60])
        g.cis = float(fileBlock[3][60:79])

        #fifth line
        g.i0 = float(fileBlock[4][0:22])
        g.crc = float(fileBlock[4][22:41])
        g.omega = float(fileBlock[4][41:60])
        g.omegaDot = float(fileBlock[4][60:79])

        #sixth line
        g.idot = float(fileBlock[5][0:22])
        g.codesL2Channel = float(fileBlock[5][22:41])
        g.gpsWeek = float(fileBlock[5][41:60])
        g.l2PDataFlag = float(fileBlock[5][60:79])

        #seventh line
        g.svAccuracy = float(fileBlock[6][0:22])
        g.svHealth = float(fileBlock[6][22:41])
        g.tgd = float(fileBlock[6][41:60])
        g.iodc = float(fileBlock[6][60:79])

        #eight Line
        g.transTime = float(fileBlock[7][0:22])
        #self.fitInterval = float(fileBlock[7][22:41])
        return g

    def _gpsFromRinex3(fileBlock):
        g = GpsMessageFile()
        #first line
        g.sat_number = fileBlock[0][0:3]
        g.epochYear = fileBlock[0][4:8]
        g.epochMonth = fileBlock[0][9:11]
        g.epochDay = fileBlock[0][12:14]
        g.epochHour = fileBlock[0][15:17]
        g.epochMinute = fileBlock[0][18:20]
        g.epochSecond = fileBlock[0][21:23]
        g.sv_clock_bias = float(fileBlock[0][24:42])
        g.sv_clock_drift = float(fileBlock[0][43:61])
        g.sv_clock_drift_rate = float(fileBlock[0][62:80])

        #second line
        g.iode = float(fileBlock[1][0:23])
        g.crs = float(fileBlock[1][23:42])
        g.deltaN = float(fileBlock[1][42:61])
        g.m0 = float(fileBlock[1][61:80])

        #third line
        g.cuc = float(fileBlock[2][0:23])
        g.eccentricity = float(fileBlock[2][23:42])
        g.cus = float(fileBlock[2][42:61])
        g.sqrtA = float(fileBlock[2][61:80])

        #fourth line
        g.toeTime = float(fileBlock[3][0:23])
        g.cic = float(fileBlock[3][23:42])
        g.omega0 = float(fileBlock[3][42:61])
        g.cis = float(fileBlock[3][61:80])

        #fifth line
        g.i0 = float(fileBlock[4][0:23])
        g.crc = float(fileBlock[4][23:42])
        g.omega = float(fileBlock[4][42:61])
        g.omegaDot = float(fileBlock[4][61:80])

        #sixth line
        g.idot = float(fileBlock[5][0:23])
        g.codesL2Channel = float(fileBlock[5][23:42])
        g.gpsWeek = float(fileBlock[5][42:61])
        g.l2PDataFlag = float(fileBlock[5][61:80])

        #seventh line
        g.svAccuracy = float(fileBlock[6][0:23])
        g.svHealth = float(fileBlock[6][23:42])
        g.tgd = float(fileBlock[6][42:61])
        g.iodc = float(fileBlock[6][61:80])

        #eight Line
        g.transTime = float(fileBlock[7][0:23])
        g.fitInterval = float(fileBlock[7][23:42])

        return g
