import math
import numpy as np
from numpy.linalg import inv
from numpy.linalg import norm

from copy import deepcopy
from copy import copy
class RinexFileReader:
    def transformGPSListDictionary(self,lista):
        dic = {}
        for gps in lista:
            if gps.sat_number in dic:
                dic[gps.sat_number].append(gps)
            else:
                dic[gps.sat_number] = [gps]
        return dic

    def readReceptorFromObservation(self,fileName):
        l = []
        with open(fileName,'r') as file:
            #print('Abriu Arquivo')
            content = file.readlines()
            file.seek(0)
            endOfHeaderLineNumber = 0
            _version = 0
            _aprox_position = (0,0,0)
            _sat_pseudoDistance = {}
            for i in range(len(content)):
                #print('Linha ',str(i),content[i])
                #Detecta a versão do arquivo rinex, atualmente aceita versão 2 ou 3
                if content[i][-21:].strip() == 'RINEX VERSION / TYPE':
                    _version = int(round(float(content[i][:10].strip()),0))
                    #print('detectou versao',str(_version))
                #Detecta o fim do header do arquivo para possibilitar a importação dos dados
                #dos GPSs
                if content[i][-21:].strip() == 'APPROX POSITION XYZ':
                    dist = content[i][:50].strip().split(' ')
                    _aprox_position = (float(dist[0]),float(dist[1]),float(dist[2]))
                if content[i][-21:].strip() == 'PRN / # OF OBS':
                    _sat_pseudoDistance[content[i][4:6].strip()] = int(content[i][7:13].strip())
                if 'END OF HEADER' in content[i]:
                    endOfHeaderLineNumber = i
                    #print('detectou header',str(i))
                    break
            if _version == 2:

                #print(content[endOfHeaderLineNumber + 1])
                total_lines = len(content) - endOfHeaderLineNumber - 1
                lines_visited = 0
                while lines_visited < total_lines:
                    #print(lines_visited,total_lines)
                    linhas_pular = int(content[endOfHeaderLineNumber + 1 + lines_visited][30:32])
                    #print(linhas_pular,content[endOfHeaderLineNumber + 1 + lines_visited])

                    gpsData = content[endOfHeaderLineNumber + 1 + lines_visited:endOfHeaderLineNumber  + lines_visited+linhas_pular+2]
                    lines_visited += linhas_pular+1
                    for i in range(len(gpsData)):
                        gpsData[i] = gpsData[i].replace('D', 'E')
                    if _version == 2:
                        l.append(GPSFactory.createReceptorFromRinexFile(gpsData, _aprox_position, _sat_pseudoDistance, _version))
            else:
                l = None
        return l
    def readGPSFromEphemeris(self,fileName):
        l = []

        with open(fileName,'r') as file:
            #print('Abriu Arquivo')
            content = file.readlines()
            file.seek(0)
            endOfHeaderLineNumber = 0
            _version = 0
            for i in range(len(content)):
                #print('Linha ',str(i),content[i])
                #Detecta a versão do arquivo rinex, atualmente aceita versão 2 ou 3
                if content[i][-21:].strip() == 'RINEX VERSION / TYPE':
                    _version = int(round(float(content[i][:10].strip()),0))
                    #print('detectou versao',str(_version))
                #Detecta o fim do header do arquivo para possibilitar a importação dos dados
                #dos GPSs
                if 'END OF HEADER' in content[i]:
                    endOfHeaderLineNumber = i
                    #print('detectou header',str(i))
                    break
            if _version == 2 or _version == 3:
                iterationNumber = (len(content) - endOfHeaderLineNumber-1) // 8
                for i in range(iterationNumber):

                    posToCut = endOfHeaderLineNumber+1+(8*i)
                    gpsData = content[posToCut:posToCut+9]
                    for i in range(len(gpsData)):
                        gpsData[i] = gpsData[i].replace('D', 'E')
                    if _version == 2 or gpsData[0][0] == 'G':
                        l.append(GPSFactory.createGPSFromRinexFile(gpsData,_version))
            else:
                l = None
        return l

class GpsMessageFile:
    def __init__(self):
        """diffTimeEphObs e o tempo que passou desde o horario importado que consta no arquivo importado"""
        self.diffTimeEphObs = 0

        """PD e a informacao de pseudo distancia que deve ser passada para tornar os calculos de posicao mais precisos"""
        self.PD = 0
        """o atributo coord_mult serve para aplicar um multiplicador 10**coord_mult ao resultado, permitindo variar
        a unidade da resposta entre metros e quilometros, por exemplo"""
        self.coord_mult = 0
    def deltaTsv(self):
        #PRECISA DE REVISAO
        #%tr   = (dia*24+hora)*3600+(minu*60)+seg;
        #tr=(0*24+0)*3600+(0*60)+30
        tr = self.toeTime
        #PD = 0
        #diffTimeEphObs = 0
        c = 299792458
        tgps = tr - self.PD/c;

        dts  = self.sv_clock_bias + self.sv_clock_drift*(tgps - self.toeTime) + self.sv_clock_drift_rate*(tgps-self.toeTime)**2;

        #%tempo de propagação aproximado do sinal
        tal = self.PD/c - self.diffTimeEphObs + dts

        tgps = tr - tal + dts;

        #%time diff from toe

        dt = tgps - self.toeTime;
        return dt

    def F(self):
        return -4.442807633E-10 # Constant, [sec/(meter)^(1/2)]
    def GM(self):
        return 3.986005E+14
    def omegaEarth(self):
        return 7.292115090E-5
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
        while e0 != e1 and lim < 15:
            e0 = e1
            e1 = self.M() + self.eccentricity*math.sin(e0)
            lim += 1
            #print('lim',lim,e1)
        return e1
    """dtr relativistc correction"""
    def dtr(self):
        return self.F() * self.eccentricity * self.sqrtA * math.sin(self.E())
    def Rel(self):
        return self.dtr()*299792458
    """Compute satellite clock correction"""
    def julianDay(self):
        # Converts the two digit year to a four digit year.
        # Two digit year represents a year in the range 1980-2079.
        if self.epochYear >= 80 and self.epochYear <= 99:
            year = 1900 + self.epochYear;

        if self.epochYear >= 0 and self.epochYear <= 79:
            year = 2000 + self.epochYear;

        # Calculates the 'm' term used below from the given calendar month.
        if self.epochMonth <= 2:
            y = year - 1;
            m = self.epochMonth + 12;
        if self.epochMonth > 2:
            y = year
            m = self.epochMonth

        #Computes the Julian date corresponding to the given calendar date.
        return  round((365.25 * y),0) + round((30.6001 * (m+1)),0) + self.epochDay + ( (self.epochHour + self.epochMinute / 60 + self.epochSecond / 3600) / 24 ) + 1720981.5
    def gpsWeek2(self):
        return round( (self.julianDay() - 2444244.5) / 7 ,0)
    def gpsSeconds(self):
        secs_per_week = 604800
        return round(((((self.julianDay()-2444244.5) / 7)-self.gpsWeek2())*secs_per_week)/0.5,0)*0.5
    def satClkCorr(self):



        clkCorr = (self.sv_clock_drift_rate * (self.diffTimeEphObs) + self.sv_clock_drift) * (self.diffTimeEphObs) + self.sv_clock_bias

        return clkCorr*299792458
    #True Anomaly
    def V(self):
        senv = (((1-self.eccentricity**2)**0.5)*math.sin(self.E()))/(1-self.eccentricity*math.cos(self.E()))
        cosv = (math.cos(self.E())-self.eccentricity)/(1-self.eccentricity*math.cos(self.E()))
        tanv = senv/cosv

        #atanv = math.atan(tanv);
        #if senv >= 0.0 and cosv <= 0.0:
        #    atanv = math.pi + math.atan(tanv)
        #elif senv < 0.0 and cosv < 0.0:
        #    atanv = math.pi + math.atan(tanv)
        #elif senv < 0.0 and cosv >= 0.0:
        #    atanv = 2*math.pi + math.atan(tanv);

        return math.atan2(senv,cosv)

        #print('eccentricity',str(self.eccentricity))
        #print('E',str(self.E()))
        #print(str(math.acos((math.cos(self.E())-self.eccentricity)/(1-self.eccentricity*math.cos(self.E())))))
        #print(str(math.asin((((1-self.eccentricity**2)**0.5)*math.sin(self.E()))/(1-self.eccentricity*math.cos(self.E())))))
        #print(self.sat_number,str(senv),str(cosv),str(tanv),str((math.atan((((1-self.eccentricity**2)**0.5)*math.sin(self.E()))/(math.cos(self.E())-self.eccentricity)))))
        return  (math.atan((((1-self.eccentricity**2)**0.5)*math.sin(self.E()))/(math.cos(self.E())-self.eccentricity)))


        #return math.acos((math.cos(self.E())-self.eccentricity)/(1-self.eccentricity*math.cos(self.E())))
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

class GPSReceptor:
    def __init__(self):
        self.aprox_pos = (0,0,0)
        self.sat_number = []
        self.gps = []
        self.epochYear = ''
        self.epochMonth = ''
        self.epochDay = ''
        self.epochHour = ''
        self.epochMinute = ''
        self.epochSecond = ''
        self.sat_pseudo = {}
    def epoch(self):
        return str(self.epochYear) +'-'+str(self.epochMonth)+'-'+str(self.epochDay)+'-'+str(self.epochHour)+'-'+str(self.epochMinute)+'-'+str(self.epochSecond)
    def loadSateliteData(self,dic):
        for sat in self.sat_number:
            #print('sat',sat)
            lista = dic[sat]
            if len(lista) > 0:
                lista = sorted(lista,key=lambda gps: abs(self.epochYear - gps.epochYear)*365*24*60*60 + abs(self.epochMonth - gps.epochMonth)*30*24*60*60 + abs(self.epochDay - gps.epochDay)*24*60*60 + abs(self.epochHour - gps.epochHour)*60*60 + abs(self.epochMinute - gps.epochMinute)*60 + abs(self.epochSecond - gps.epochSecond))
                g = copy(lista[0])
                g.CorrP1 = self.sat_obs_data[g.sat_number]['P3'] + g.satClkCorr() + g.Rel()
                g.diffTimeEphObs = (self.epochYear - g.epochYear)*365*24*60*60 + (self.epochMonth - g.epochMonth)*30*24*60*60 + (self.epochDay - g.epochDay)*24*60*60 + (self.epochHour - g.epochHour)*60*60 + (self.epochMinute - g.epochMinute)*60 + (self.epochSecond - g.epochSecond)
                #print('diffTimeEphObs',g.diffTimeEphObs)
                g.PD = self.sat_pseudo[g.sat_number]
                self.gps.append(g)
                #pseudodistancia
    def getCoordinates(self):
        if len(self.gps) <= 0:
            return None
        lista_spos = []
        lista_a = []
        lista_b =  []
        corrP1 = []
        for gps in self.gps:

            #print('gps',gps.sat_number,str(gps.epochHour),str(gps.epochMinute),str(gps.epochSecond),'diffTimeEphObs',gps.diffTimeEphObs)
            gps_coord = gps.coordinate_WGS84()
            lista_spos.append([gps_coord[0],gps_coord[1],gps_coord[2]])
            corrP1.append(gps.CorrP1)
            lista_b.append(0)

        matrix_b = np.matrix(lista_b)

        matrix_spos = np.matrix(lista_spos)

        rec_xyz = np.matrix([self.aprox_pos[0],self.aprox_pos[1],self.aprox_pos[2]])
        #for j in range(10):
        for i in range(len(self.gps)):

            matrix_b[0,i] = (corrP1[i] - norm(matrix_spos[i] - rec_xyz[0] ,ord='fro'))
            lista_a.append([(-(matrix_spos[i,0] - rec_xyz[0,0])) / corrP1[i] ,
                            (-(matrix_spos[i,1] - rec_xyz[0,1])) / corrP1[i] ,
                            (-(matrix_spos[i,2] - rec_xyz[0,2])) / corrP1[i] ,
                            1 ])
        matrix_a = np.matrix(lista_a)
        delta_xyzb = np.linalg.lstsq(matrix_a,matrix_b.transpose())[0]
        rec_xyz = rec_xyz - delta_xyzb[:-1].transpose()
        lista_a = []
        return  (rec_xyz[0,0],rec_xyz[0,1],rec_xyz[0,2])
    def getCoordinates_(self):
        if len(self.gps) <= 0:
            return None

        lista_X0 = [self.aprox_pos[0],self.aprox_pos[1],self.aprox_pos[2],0]
        #lista_X0 = [0,0,0,0]
        matrix_X0 = np.matrix(lista_X0).transpose()

        for i in range(1):
            lista_a = []
            lista_L = []

            for gps in self.gps:

                #print('gps',gps.sat_number,str(gps.epochHour),str(gps.epochMinute),str(gps.epochSecond),'diffTimeEphObs',gps.diffTimeEphObs)
                gps_coord = gps.coordinate_WGS84()
                #print(gps_coord)
                denominator = ((lista_X0[0] - gps_coord[0])**2 + (lista_X0[1] - gps_coord[1])**2 + (lista_X0[2] - gps_coord[2])**2)**0.5
                lista_a.append([(lista_X0[0] - gps_coord[0])/denominator,
                                (lista_X0[1] - gps_coord[1])/denominator,
                                (lista_X0[2] - gps_coord[2])/denominator,
                                1])
                lista_L.append(self.sat_pseudo[gps.sat_number])
            matrix_L = np.matrix(lista_L).transpose()
            matriz_a = np.matrix(lista_a)
            matriz_aT = matriz_a.transpose()
            matrix_deltaX = np.dot(inv(np.dot(matriz_aT,matriz_a)),
                                    np.dot(matriz_aT,matrix_L))

            #matrix_deltaX[0][0] = 90
            matrix_Xa = matrix_X0 + matrix_deltaX

            print('mariz deltaX',matrix_deltaX)
            print('mariz X0',matrix_X0)
            print('mariz Xa',matrix_Xa)

            lista_X0 = [matrix_Xa.A[0][0],matrix_Xa.A[1][0],matrix_Xa.A[2][0],matrix_Xa.A[3][0]]
            matrix_X0 = np.matrix(lista_X0).transpose()
        return (matrix_Xa.A[0][0],matrix_Xa.A[1][0],matrix_Xa.A[2][0])


"""Classe responsavel por criar instancias de GpsMessageFile de acordo com o Layout Rinex 2 ou Rinex 3"""
class GPSFactory:
    def createReceptorFromRinexFile(fileBlock, aprox_pos, sat_pseudo, version):
        if version == 2:
            g = GPSFactory._receptorFromRinex2(fileBlock, aprox_pos, sat_pseudo)
        return g
    def _receptorFromRinex2(fileBlock,aprox_pos, sat_pseudo):

        fL1 = 1575.42E6   # L1 frequency (Hz)
        fL2 = 1227.6E6    # L2 frequency (Hz)
        B=fL2**2/(fL2**2-fL1**2)
        A=-B+1;


        r = GPSReceptor()
        r.aprox_pos = aprox_pos
        r.sat_number = []
        r.gps = []
        r.sat_pseudo = sat_pseudo
        r.sat_obs_data = {}
        auxString = fileBlock[0][31:].strip()
        auxSats = auxString.split('G')
        for i in range(1,len(auxSats)):
            r.sat_number.append(auxSats[i].strip())
            c1 = float(fileBlock[i][2:14].strip())
            l1 = float(fileBlock[i][17:32].strip())
            p2 = float(fileBlock[i][34:46].strip())
            l2 = float(fileBlock[i][50:].strip())

            p3=A*c1+B*p2


            r.sat_obs_data[auxSats[i].strip()] = {'C1':c1,'L1':l1,'P2':p2,'L2':l2,'P3':p3}
        print(r.sat_obs_data)
        r.epochYear = int(fileBlock[0][0:3])
        r.epochMonth = int(fileBlock[0][4:6])
        r.epochDay = int(fileBlock[0][8:10])
        r.epochHour = int(fileBlock[0][10:12])
        r.epochMinute = int(fileBlock[0][13:15])
        r.epochSecond = int(round(float(fileBlock[0][16:26]),0))
        return r

    def createGPSFromRinexFile(fileBlock, version):
        if version == 2:
            g = GPSFactory._gpsFromRinex2(fileBlock)
        else:
            g = GPSFactory._gpsFromRinex3(fileBlock)
        return g
    def _gpsFromRinex2(fileBlock):
        g = GpsMessageFile()
        #first line
        g.sat_number = fileBlock[0][0:3].strip()
        g.epochYear = int(fileBlock[0][3:5])
        g.epochMonth = int(fileBlock[0][6:8])
        g.epochDay = int(fileBlock[0][9:11])
        g.epochHour = int(fileBlock[0][12:14])
        g.epochMinute = int(fileBlock[0][15:17])
        g.epochSecond = int(round(float(fileBlock[0][18:22]),0))
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
        g.sat_number = fileBlock[0][0:3].strip()
        g.epochYear = int(fileBlock[0][4:8])
        g.epochMonth = int(fileBlock[0][9:11])
        g.epochDay = int(fileBlock[0][12:14])
        g.epochHour = int(fileBlock[0][15:17])
        g.epochMinute = int(fileBlock[0][18:20])
        g.epochSecond = int(fileBlock[0][21:23])
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
