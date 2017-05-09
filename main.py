from flask import Flask, render_template, request
from werkzeug.utils import secure_filename
import os
import sys
sys.path.append("..")
from gpsMessageFile import RinexFileReader

import folium
import pyproj

UPLOAD_FOLDER = '.'
r = RinexFileReader()
app=Flask(__name__)
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/map',methods = ['GET'])
def show_map():
    return render_template('mapa.html')

@app.route('/files',methods = ['POST'])
def files():
    lista = []
    receptores = []
    recep_comGPS = []
    max_obs =  int(request.form['maxobs'])

    #print('Caiu no servidor')
    if request.method == 'POST':
        #print('Metodo certo')
        f = request.files['nav_file']

        if f.filename == '':
            pass
        if f:
            filename = secure_filename(f.filename)
            f.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))
        #print(filename)
        lista = r.readFromEphemeris(filename)
        os.remove(os.path.join(app.config['UPLOAD_FOLDER'], filename))

        f2 = request.files['obs_file']

        if f2.filename == '':
            pass
        if f2:
            filename = secure_filename(f2.filename)
            f2.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))
        #print(filename)
            recep_comGPS = []
            receptor = r.readFromObservation(filename,lista,max_obs)
            receptor.calculate_receiverPosition()
            os.remove(os.path.join(app.config['UPLOAD_FOLDER'], filename))


            #salvaMapa(converteLatLong(receptor))
            listaLatLong = converteLatLong(receptor)
            #print(listaLatLong)
            listaGPS = []
            for sat in lista:
                for ephem in lista[sat]:
                    listaGPS.append(ephem)
            listaGPS = sorted(listaGPS,key=lambda e: e.toc())
            #print('Verificar',receptores[55].aprox_pos, receptores[55].sat_number,receptores[55].gps,receptores[55].sat_pseudo)
            #print(receptores[55].getCoordinates())
            #print('Verificou')
        #print(lista[0].coordinate_WGS84())
    return render_template('index.html',listaGPS = listaGPS, receptor = receptor,listaLatLong = listaLatLong)

def converteLatLong(receptor):
    listaLatLong = []
    for obs in receptor.observations:
        x,y,z = obs.rec_xyz
        ecef = pyproj.Proj(proj='geocent', ellps='WGS84', datum='WGS84')
        lla = pyproj.Proj(proj='latlong', ellps='WGS84', datum='WGS84')
        lon, lat, alt = pyproj.transform(ecef, lla, x, y, z, radians=False)
        listaLatLong.append((str(obs.gpsSeconds()),lat,lon,alt))
    return listaLatLong

if __name__=='__main__':
    app.run(debug=True)
