import math
import numpy as np
from numpy.linalg import inv
from numpy.linalg import norm

from copy import deepcopy
from copy import copy
from math import floor
from datetime import datetime

#Constantes uteis para os calculos
c = 299792458          #----> Speed of light (meters/s).
L1 = 1575.42E6         #----> Freqs in Hz.
L2 = 1227.60E6
L5 = 1176.45E6

LAMBDA_L1 = c / L1     #----> Wavelengths in meters.
LAMBDA_L2 = c / L2
LAMBDA_L5 = c / L5


class RinexFileReader:
    def readFromObservation(self,fileName,ephem,max_obs):
        with open(fileName,'r') as file:
            #print('Abriu Arquivo')
            content = file.readlines()
        lista_obs = []
        _aprox_position = (0,0,0)
        for i in range(len(content)):
            #Detecta a versão do arquivo rinex, atualmente aceita versão 2 ou 3
            if content[i][-21:].strip() == 'RINEX VERSION / TYPE':
                _version = int(round(float(content[i][:10].strip()),0))
                #print('detectou versao',str(_version))

            if content[i][-21:].strip() == 'APPROX POSITION XYZ':
                dist = content[i][:50].strip().split(' ')
                _aprox_position = (float(dist[0]),float(dist[1]),float(dist[2]))
            #Detecta o fim do header do arquivo para possibilitar a importação dos dados
            #dos GPSs
            if 'END OF HEADER' in content[i]:
                endOfHeaderLineNumber = i
                #print('detectou header',str(i))
                break

        if _version == 2:
            total_lines = len(content) - endOfHeaderLineNumber - 1
            lines_visited = 0
            while lines_visited < total_lines:
                linhas_pular = int(content[endOfHeaderLineNumber + 1 + lines_visited][30:32])
                fileBlock = content[endOfHeaderLineNumber + 1 + lines_visited:endOfHeaderLineNumber  + lines_visited+linhas_pular+2]
                lines_visited += linhas_pular+1
                for i in range(len(fileBlock)):
                    fileBlock[i] = fileBlock[i].replace('D', 'E')
                if _version == 2:
                    o = GPSFactory.createObservationFromRinexFile(fileBlock,_version)
                    lista_obs.append(o)
                    if max_obs == 1:
                        break
                    elif max_obs != 0:
                        max_obs -= 1

        else:
            lista_obs = None
        r = Receiver()
        r.observations = lista_obs
        r._aprox_position = _aprox_position
        r.ephemeris = ephem
        return r

    def readFromEphemeris(self,fileName):


        with open(fileName,'r') as file:
            #print('Abriu Arquivo')
            content = file.readlines()

        endOfHeaderLineNumber = 0
        _version = 0
        lista_ephe = {}
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
                fileBlock = content[posToCut:posToCut+9]
                for i in range(len(fileBlock)):
                    fileBlock[i] = fileBlock[i].replace('D', 'E')
                if _version == 2 or fileBlock[0][0] == 'G':
                    g = GPSFactory.createEphemerisFromRinexFile(fileBlock,_version)
                    if g.sat_number not in lista_ephe:
                        lista_ephe[g.sat_number] = []
                    lista_ephe[g.sat_number].append(g)
                else:
                    lista_ephe = None
        return lista_ephe
class SatOrbit:
    def __init__(self):
        self.sat_number = 0
        self.epochYear = 0
        self.epochMonth = 0
        self.epochDay = 0
        self.epochHour = 0
        self.epochMinute = 0
        self.epochSecond = 0
        self.c1 = 0
        self.l1 = 0
        self.p2 = 0
        self.l2 = 0
class GPSUtils:
    def gpsWeek_seconds(epochYear,epochMonth,epochDay,epochHour,epochMinute,epochSecond):
        secs_per_week = 604800
        year = epochYear
        # Converts the two digit year to a four digit year.
        # Two digit year represents a year in the range 1980-2079.
        if year >= 80 and year <= 99:
            year = 1900 + year;

        if year >= 0 and year <= 79:
            year = 2000 + year;

        # Calculates the 'm' term used below from the given calendar month.
        if epochMonth <= 2:
            y = year - 1;
            m = epochMonth + 12;
        if epochMonth > 2:
            y = year
            m = epochMonth

        #Computes the Julian date corresponding to the given calendar date.
        julianDay = floor((365.25 * y)) + floor((30.6001 * (m+1))) + epochDay + ( (epochHour + epochMinute / 60 + epochSecond / 3600) / 24 ) + 1720981.5
        gpsWeek = floor( (julianDay - 2444244.5) / 7)
        gpsSeconds = round(((((julianDay-2444244.5) / 7)-gpsWeek)*secs_per_week)/0.5,0)*0.5
        return (gpsWeek,gpsSeconds)
    def kepOrb2E(M,e):
        #print('Começo Kepler')
        if -math.pi < M < 0 or M > math.pi:
            E = M - e;
        else:
            E = M + e;


        check = 1
        i = 10
        while check > 10E-10 and i > 0:
            #print('teste',check,i)
            E_new = (E + (M - E + e * math.sin(E))/(1 - e * math.cos(E)))
            check = abs(E_new - E)
            E = E_new
            i -= 1
        #print('Fim Kepler')
        return E

class Receiver:
    def __init__(self):
        self.observations = []
        self.ephemeris = {}
        self._aprox_position = (0,0,0)

    def get_broadcast_orbits(self,obs,rec_pos):
        meu = 3.986005E14;         # earth's universal gravitational [m^3/s^2]
        odote = 7.2921151467E-5;   # earth's rotation rate (rad/sec)
        lightspeed = 2.99792458E8; # speed of light (m/s)

        F = -4.442807633E-10; # Constant, [sec/(meter)^(1/2)]

        tsat = obs.gpsSeconds()
        #print('Sat #',obs.sat_number)
        for j in range(4):
            lista = self.ephemeris[obs.sat_number]
            lista = sorted(lista,key=lambda e: abs(tsat - e.gpsWeek()))
            ephem = copy(lista[0])
            ephem.a = ephem.sqrtA **2
            ephem.dn = ephem.delta_n
            ephem.omg0 = ephem.OMEGA
            ephem.odot = ephem.OMEGA_dot

            n0 = (meu/ephem.a**3)**0.5
            #print('n0',n0)
            t = tsat-ephem.toe
            n = n0 + ephem.dn
            m = ephem.M0 + n*t
            #print('M0 coisas',ephem.M0,n,t,tsat,ephem.toe)
            m_dot=n #% Calculate Velocity
            E = GPSUtils.kepOrb2E(m,ephem.e)
            #print('E',E)

            #Compute relativistic correction term
            dtr = F * ephem.e * ((ephem.a)**0.5) * math.sin(E)
            obs.Rel = dtr*c

            # Compute satellite clock correction
            clkCorr= (ephem.af2 * (tsat-ephem.toc()) + ephem.af1) * (tsat-ephem.toc()) + ephem.af0
            obs.satClkCorr = clkCorr*c
            #print('ClockCorr',obs.satClkCorr)
            t = t - clkCorr;

            E_dot=m_dot/(1-ephem.e*math.cos(E)) #Calculate Velocity

            obs.E = E
            v = math.atan2(((1-ephem.e**2)**0.5)*math.sin(E),math.cos(E)-ephem.e)

            v_dot=math.sin(E)*E_dot*(1+ephem.e*math.cos(v))/(math.sin(v)*(1-ephem.e*math.cos(E))) # Calculate Velocity


            phi = v + ephem.omega

            phi_dot=v_dot                            # Calculate Velocity

            du = ephem.Cus*math.sin(2*phi) + ephem.cuc*math.cos(2*phi)
            dr = ephem.crs*math.sin(2*phi) + ephem.crc*math.cos(2*phi)
            di = ephem.Cis*math.sin(2*phi) + ephem.Cic*math.cos(2*phi)

            du_dot=2*(ephem.Cus*math.cos(2*phi)-ephem.cuc*math.sin(2*phi))*phi_dot # Calculate Velocity
            dr_dot=2*(ephem.crs*math.cos(2*phi)-ephem.crc*math.sin(2*phi))*phi_dot # Calculate Velocity
            di_dot=2*(ephem.Cis*math.cos(2*phi)-ephem.Cic*math.sin(2*phi))*phi_dot # Calculate Velocity



            u = phi + du
            r = ephem.a*(1-ephem.e*math.cos(E)) + dr
            i = ephem.i0 + di + ephem.i_dot*t

            u_dot=phi_dot+du_dot#              % Calculate Velocity
            r_dot=ephem.a*ephem.e*math.sin(E)*E_dot+dr_dot#  % Calculate Velocity
            i_dot=ephem.i_dot+di_dot#                     % Calculate Velocity

            xp = r*math.cos(u)
            yp = r*math.sin(u)

            xp_dot=r_dot*math.cos(u)-r*math.sin(u)*u_dot#          % Calculate Velocity
            yp_dot=r_dot*math.sin(u)+r*math.cos(u)*u_dot#           % Calculate Velocity

            omg = ephem.omg0 + (ephem.odot - odote)*t - odote*ephem.toe

            omg_dot=ephem.OMEGA_dot - odote                 # % Calculate Velocity

            obs.XS = xp*math.cos(omg) - yp*math.cos(i)*math.sin(omg)
            obs.YS = xp*math.sin(omg) + yp*math.cos(i)*math.cos(omg)
            obs.ZS = yp*math.sin(i)

            obs.VXS= xp_dot*math.cos(omg)-yp_dot*math.cos(i)*math.sin(omg)+yp*math.sin(i)*math.sin(omg)*i_dot-obs.YS*omg_dot
            obs.VYS= xp_dot*math.sin(omg)+yp_dot*math.cos(i)*math.cos(omg)-yp*math.sin(i)*i_dot*math.cos(omg)+obs.XS*omg_dot
            obs.VZS= yp_dot*math.sin(i)+yp*math.cos(i)*i_dot

            # compute the range
            R = ((obs.ZS - rec_pos[0])**2 + (obs.ZS - rec_pos[1])**2 + (obs.ZS - rec_pos[2])**2)**0.5

            obs.RANGE = R    #% range
            tau=R/lightspeed;

            #% Add earth rotation correction here
            phi=-odote*tau;
            l1 = np.matrix([[math.cos(phi),-math.sin(phi)],[math.sin(phi),math.cos(phi)]])
            l2 = np.matrix([obs.XS,obs.YS])

            corr= np.dot(l1,l2.T)
            #print(corr[0,0])
            obs.XS=corr[0,0]
            obs.YS=corr[1,0]
            #% update the range
            R_new = ((obs.XS - rec_pos[0])**2 + (obs.YS - rec_pos[1])**2 + (obs.ZS - rec_pos[2])**2)**0.5
            obs.RANGE = R_new    #% range
            tau_new=R_new/lightspeed;

            tsat = obs.gpsSeconds() - tau_new;
            #print('tsat apos',tsat,tau_new)
    def calculate_receiverPosition(self):
        fL1 = 1575.42E6   # L1 frequency (Hz)
        fL2 = 1227.6E6    # L2 frequency (Hz)
        B=fL2**2/(fL2**2-fL1**2)
        A=-B+1
        for o in self.observations:
            rec_xyz = self._aprox_position
            for obs_data in o.data:
                sat = SatOrbit()
                sat.sat_number = obs_data.sat_number
                sat.c1 = obs_data.c1
                sat.l1 = obs_data.l1
                sat.p2 = obs_data.p2
                sat.l2 = obs_data.l2
                sat.p3 = A*sat.c1+B*sat.p2
                obs_data.sat = sat

            for i in range(10):
                for obs_data in o.data:
                    self.get_broadcast_orbits(obs_data,rec_xyz)
                    obs_data.clk = obs_data.satClkCorr
                    obs_data.CorrP1 = obs_data.sat.p3 + obs_data.clk + obs_data.Rel
                    #print('Sat! ',obs_data.sat_number)
                    #print('Corrp1',obs_data.CorrP1)
                    #print('P3',obs_data.sat.p3)
                    #print('clk',obs_data.clk)
                    #print('Rel',obs_data.Rel)


                dxyz = self.delta_xyz(o,rec_xyz)
                #print(o.gpsSeconds(),i,dxyz)
                rec_xyz = rec_xyz + dxyz
                rec_xyz = (rec_xyz[0,0],rec_xyz[0,1],rec_xyz[0,2])
                #print('REC_XYZ',rec_xyz)
                o.rec_xyz = rec_xyz
    def delta_xyz(self,observation,aprox_pos):
        lista_spos = []
        lista_a = []
        lista_b =  []
        corrP1 = []
        for obs_data in observation.data:
            lista_spos.append([obs_data.XS,obs_data.YS,obs_data.ZS])
            corrP1.append(obs_data.CorrP1)
            lista_b.append(0)

        matrix_b = np.matrix(lista_b)

        matrix_spos = np.matrix(lista_spos)

        rec_xyz = np.matrix([aprox_pos[0],aprox_pos[1],aprox_pos[2]])
        #for j in range(10):
        for i in range(len(observation.data)):

            matrix_b[0,i] = (corrP1[i] - norm(matrix_spos[i] - rec_xyz ,ord='fro'))
            #print('Corrp1',corrP1[i])
            #print('matrix SPOPS',matrix_spos[i])
            #print('recxyz',rec_xyz)

            #print(matrix_b)
            lista_a.append([(-(matrix_spos[i,0] - rec_xyz[0,0])) / corrP1[i] ,
                            (-(matrix_spos[i,1] - rec_xyz[0,1])) / corrP1[i] ,
                            (-(matrix_spos[i,2] - rec_xyz[0,2])) / corrP1[i] ,
                            1 ])
        matrix_a = np.matrix(lista_a)
        delta_xyzb = np.linalg.lstsq(matrix_a,matrix_b.transpose())[0]
        return delta_xyzb[:-1].transpose()



class Ephemeris:
    def _gpsWeekSeconds(self):
        return GPSUtils.gpsWeek_seconds(self.epochYear,self.epochMonth,self.epochDay,
                                        self.epochHour,self.epochMinute,self.epochSecond)
    def gpsWeek(self):
        return self._gpsWeekSeconds()[0]
    def toc(self):
        return self._gpsWeekSeconds()[1]
    def calculatePosition(self):
        meu = 3.986005E14;         # earth's universal gravitational [m^3/s^2]
        odote = 7.2921151467E-5;   # earth's rotation rate (rad/sec)
        lightspeed = 2.99792458E8; # speed of light (m/s)

        F = -4.442807633E-10; # Constant, [sec/(meter)^(1/2)]

        tsat = self.toc()
        #print('Sat #',obs.sat_number)



        self.a = self.sqrtA **2
        self.dn = self.delta_n
        self.omg0 = self.OMEGA
        self.odot = self.OMEGA_dot

        n0 = (meu/self.a**3)**0.5
        #print('n0',n0)
        t = tsat-self.toe
        n = n0 + self.dn
        m = self.M0 + n*t
        #print('M0 coisas',self.M0,n,t,tsat,self.toe)
        m_dot=n #% Calculate Velocity
        E = GPSUtils.kepOrb2E(m,self.e)
        #print('E',E)

        #Compute relativistic correction term
        dtr = F * self.e * ((self.a)**0.5) * math.sin(E)
        Rel = dtr*c

        # Compute satellite clock correction
        clkCorr= (self.af2 * (tsat-self.toc()) + self.af1) * (tsat-self.toc()) + self.af0
        satClkCorr = clkCorr*c
        #print('ClockCorr',obs.satClkCorr)
        t = t - clkCorr;

        E_dot=m_dot/(1-self.e*math.cos(E)) #Calculate Velocity


        v = math.atan2(((1-self.e**2)**0.5)*math.sin(E),math.cos(E)-self.e)

        v_dot=math.sin(E)*E_dot*(1+self.e*math.cos(v))/(math.sin(v)*(1-self.e*math.cos(E))) # Calculate Velocity


        phi = v + self.omega

        phi_dot=v_dot                            # Calculate Velocity

        du = self.Cus*math.sin(2*phi) + self.cuc*math.cos(2*phi)
        dr = self.crs*math.sin(2*phi) + self.crc*math.cos(2*phi)
        di = self.Cis*math.sin(2*phi) + self.Cic*math.cos(2*phi)

        du_dot=2*(self.Cus*math.cos(2*phi)-self.cuc*math.sin(2*phi))*phi_dot # Calculate Velocity
        dr_dot=2*(self.crs*math.cos(2*phi)-self.crc*math.sin(2*phi))*phi_dot # Calculate Velocity
        di_dot=2*(self.Cis*math.cos(2*phi)-self.Cic*math.sin(2*phi))*phi_dot # Calculate Velocity



        u = phi + du
        r = self.a*(1-self.e*math.cos(E)) + dr
        i = self.i0 + di + self.i_dot*t

        u_dot=phi_dot+du_dot#              % Calculate Velocity
        r_dot=self.a*self.e*math.sin(E)*E_dot+dr_dot#  % Calculate Velocity
        i_dot=self.i_dot+di_dot#                     % Calculate Velocity

        xp = r*math.cos(u)
        yp = r*math.sin(u)

        xp_dot=r_dot*math.cos(u)-r*math.sin(u)*u_dot#          % Calculate Velocity
        yp_dot=r_dot*math.sin(u)+r*math.cos(u)*u_dot#           % Calculate Velocity

        omg = self.omg0 + (self.odot - odote)*t - odote*self.toe

        omg_dot=self.OMEGA_dot - odote                 # % Calculate Velocity

        XS = xp*math.cos(omg) - yp*math.cos(i)*math.sin(omg)
        YS = xp*math.sin(omg) + yp*math.cos(i)*math.cos(omg)
        ZS = yp*math.sin(i)
        return (XS,YS,ZS)
class Observation:
    def __init__(self):
        self.sat_number = 0
        self.epochYear = 0
        self.epochMonth = 0
        self.epochDay = 0
        self.epochHour = 0
        self.epochMinute = 0
        self.epochSecond = 0
        self.data = []

    def getStrDate(self):
        stringEpoch = ('00'+str(self.epochYear))[-2:] + ' '+ str(self.epochMonth) + ' ' + str(self.epochDay) + ' '+ str(self.epochHour) + ' ' + str(self.epochMinute) + ' ' + str(self.epochSecond)
        d = datetime.strptime(stringEpoch,'%y %m %d %H %M %S')
        return d.__str__()
    def _gpsWeekSeconds(self):
        return GPSUtils.gpsWeek_seconds(self.epochYear,self.epochMonth,self.epochDay,
                                        self.epochHour,self.epochMinute,self.epochSecond)
    def gpsWeek(self):
        return self._gpsWeekSeconds()[0]
    def gpsSeconds(self):
        return self._gpsWeekSeconds()[1]

class ObsData:
    def __init__(self):
        self.sat_number = 0
        self.epochYear = 0
        self.epochMonth = 0
        self.epochDay = 0
        self.epochHour = 0
        self.epochMinute = 0
        self.epochSecond = 0
        self.c1 = 0
        self.l1 = 0
        self.p2 = 0
        self.l2 = 0
    def _gpsWeekSeconds(self):
        return GPSUtils.gpsWeek_seconds(self.epochYear,self.epochMonth,self.epochDay,
                                        self.epochHour,self.epochMinute,self.epochSecond)
    def gpsWeek(self):
        return self._gpsWeekSeconds()[0]
    def gpsSeconds(self):
        return self._gpsWeekSeconds()[1]
    def __repr__(self):
        return self.__str__()
    def __str__(self):
        return str(self.sat_number) +' '+ str(self.epochYear) +' '+ str(self.epochMonth) +' '+ str(self.epochDay) +' '+ str(self.epochHour) +' '+ str(self.epochMinute) +' '+ str(self.epochSecond) +' '+ str(self.c1) +' '+ str(self.l1) +' '+ str(self.p2) +' '+ str(self.l2)












"""Classe responsavel por criar instancias de GpsMessageFile de acordo com o Layout Rinex 2 ou Rinex 3"""
class GPSFactory:
    def createObservationFromRinexFile(fileBlock, version):
        if version == 2:
            g = GPSFactory._observationFromRinex2(fileBlock)
        return g
    def _observationFromRinex2(fileBlock):
        epochYear = int(fileBlock[0][0:3])
        epochMonth = int(fileBlock[0][4:6])
        epochDay = int(fileBlock[0][8:10])
        epochHour = int(fileBlock[0][10:12])
        epochMinute = int(fileBlock[0][13:15])
        epochSecond = int(round(float(fileBlock[0][16:26]),0))

        auxString = fileBlock[0][31:].strip()
        auxSats = auxString.split('G')
        o = Observation()
        o.epochYear = epochYear
        o.epochMonth = epochMonth
        o.epochDay = epochDay
        o.epochHour = epochHour
        o.epochMinute = epochMinute
        o.epochSecond = epochSecond

        for j in range(1,len(auxSats)):
            obs = ObsData()
            obs.epochYear = epochYear
            obs.epochMonth = epochMonth
            obs.epochDay = epochDay
            obs.epochHour = epochHour
            obs.epochMinute = epochMinute
            obs.epochSecond = epochSecond


            obs.sat_number = int(auxSats[j].strip())
            obs.c1 = float(fileBlock[j][2:14].strip())
            obs.l1 = float(fileBlock[j][17:32].strip()) * LAMBDA_L1
            obs.p2 = float(fileBlock[j][34:46].strip())
            obs.l2 = float(fileBlock[j][50:].strip())
            o.data.append(obs)
        return o

    def createEphemerisFromRinexFile(fileBlock, version):
        if version == 2:
            g = GPSFactory._ephemerisFromRinex2(fileBlock)
        else:
            g = GPSFactory._ephemerisFromRinex3(fileBlock)
        return g
    def _ephemerisFromRinex2(fileBlock):
        g = Ephemeris()
        #first line
        g.sat_number = int(fileBlock[0][0:3].strip())
        g.epochYear = int(fileBlock[0][3:5])
        g.epochMonth = int(fileBlock[0][6:8])
        g.epochDay = int(fileBlock[0][9:11])
        g.epochHour = int(fileBlock[0][12:14])
        g.epochMinute = int(fileBlock[0][15:17])
        g.epochSecond = int(round(float(fileBlock[0][18:22]),0))
        g.af0 = float(fileBlock[0][22:41])
        g.af1 = float(fileBlock[0][41:60])
        g.af2 = float(fileBlock[0][60:79])

        #second line
        g.IODE = float(fileBlock[1][0:22])
        g.crs = float(fileBlock[1][22:41])
        g.delta_n = float(fileBlock[1][41:60])
        g.M0 = float(fileBlock[1][60:79])

        #third line
        g.cuc = float(fileBlock[2][0:22])
        g.e = float(fileBlock[2][22:41])
        g.Cus = float(fileBlock[2][41:60])
        g.sqrtA = float(fileBlock[2][60:79])

        #fourth line
        g.toe = float(fileBlock[3][0:22])
        g.Cic = float(fileBlock[3][22:41])
        g.OMEGA = float(fileBlock[3][41:60])
        g.Cis = float(fileBlock[3][60:79])

        #fifth line
        g.i0 = float(fileBlock[4][0:22])
        g.crc = float(fileBlock[4][22:41])
        g.omega = float(fileBlock[4][41:60])
        g.OMEGA_dot = float(fileBlock[4][60:79])

        #sixth line
        g.i_dot = float(fileBlock[5][0:22])
        g.L2_codes = float(fileBlock[5][22:41])
        g.GPS_wk = float(fileBlock[5][41:60])
        g.l2_dataflag = float(fileBlock[5][60:79])

        #seventh line
        g.SV_Acc = float(fileBlock[6][0:22])
        g.SV_health = float(fileBlock[6][22:41])
        g.TGD = float(fileBlock[6][41:60])
        g.IODC = float(fileBlock[6][60:79])

        #eight Line
        g.transTime = float(fileBlock[7][0:22])
        return g
    def _ephemerisFromRinex3(fileBlock):
        g = Ephemeris()
        #first line
        g.sat_number = fileBlock[0][0:3].strip()
        g.epochYear = int(fileBlock[0][4:8])
        g.epochMonth = int(fileBlock[0][9:11])
        g.epochDay = int(fileBlock[0][12:14])
        g.epochHour = int(fileBlock[0][15:17])
        g.epochMinute = int(fileBlock[0][18:20])
        g.epochSecond = int(fileBlock[0][21:23])
        g.af0 = float(fileBlock[0][24:42])
        g.af1 = float(fileBlock[0][43:61])
        g.af2 = float(fileBlock[0][62:80])

        #second line
        g.IODE = float(fileBlock[1][0:23])
        g.crs = float(fileBlock[1][23:42])
        g.delta_n = float(fileBlock[1][42:61])
        g.M0 = float(fileBlock[1][61:80])

        #third line
        g.cuc = float(fileBlock[2][0:23])
        g.e = float(fileBlock[2][23:42])
        g.Cus = float(fileBlock[2][42:61])
        g.sqrtA = float(fileBlock[2][61:80])

        #fourth line
        g.toe = float(fileBlock[3][0:23])
        g.Cic = float(fileBlock[3][23:42])
        g.OMEGA = float(fileBlock[3][42:61])
        g.Cis = float(fileBlock[3][61:80])

        #fifth line
        g.i0 = float(fileBlock[4][0:23])
        g.crc = float(fileBlock[4][23:42])
        g.omega = float(fileBlock[4][42:61])
        g.OMEGA_dot = float(fileBlock[4][61:80])

        #sixth line
        g.i_dot = float(fileBlock[5][0:23])
        g.L2_codes = float(fileBlock[5][23:42])
        g.GPS_wk = float(fileBlock[5][42:61])
        g.l2_dataflag = float(fileBlock[5][61:80])

        #seventh line
        g.SV_Acc = float(fileBlock[6][0:23])
        g.SV_health = float(fileBlock[6][23:42])
        g.TGD = float(fileBlock[6][42:61])
        g.IODC = float(fileBlock[6][61:80])

        #eight Line
        g.transTime = float(fileBlock[7][0:23])
        g.fitInterval = float(fileBlock[7][23:42])

        return g
