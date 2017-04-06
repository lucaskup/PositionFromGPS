from flask import Flask, render_template, request
from werkzeug.utils import secure_filename
import os
import sys
sys.path.append("..")
from gpsMessageFile import RinexFileReader


UPLOAD_FOLDER = '.'
r = RinexFileReader()
app=Flask(__name__)
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/files',methods = ['POST'])
def files():
    lista = []
    receptores = []
    recep_comGPS = []
    max_receptor_count =  int(request.form['maxobs'])

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
        lista = r.readGPSFromEphemeris(filename)
        os.remove(os.path.join(app.config['UPLOAD_FOLDER'], filename))

        f2 = request.files['obs_file']

        if f2.filename == '':
            pass
        if f2:
            filename = secure_filename(f2.filename)
            f2.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))
        #print(filename)
            recep_comGPS = []
            receptores = r.readReceptorFromObservation(filename)
            os.remove(os.path.join(app.config['UPLOAD_FOLDER'], filename))
            dic =  r.transformGPSListDictionary(lista)
            #receptores[55].loadSateliteData(dic)
            for i in range(len(receptores)):
                if max_receptor_count != None and i>= max_receptor_count:
                    break
                receptores[i].loadSateliteData(dic)
                recep_comGPS.append(receptores[i])

            #print('Verificar',receptores[55].aprox_pos, receptores[55].sat_number,receptores[55].gps,receptores[55].sat_pseudo)
            #print(receptores[55].getCoordinates())
            #print('Verificou')
        #print(lista[0].coordinate_WGS84())
    return render_template('index.html',listaGPS = lista, listaReceptores = recep_comGPS)


if __name__=='__main__':
    app.run(debug=True)
