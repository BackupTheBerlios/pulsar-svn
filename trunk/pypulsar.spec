__version__='0.1.1'
__console__=True
installpath='D:\\program_files\\Pulsar\\trunk\\'
imgpath=installpath+'images\\'

a = Analysis(
	[os.path.join(HOMEPATH,'support\\_mountzlib.py'), 
	os.path.join(HOMEPATH,'support\\useUnicode.py'), 
	installpath+'pypulsar.py'],     
	pathex=['D:\\program_files\\Pulsar\\pyinstaller'])

pyz = PYZ(a.pure)

exe = EXE(pyz,
          a.scripts,
          exclude_binaries=1,
          name='buildpypulsar\\pypulsar.exe',
          debug=False,
          strip=False,
          upx=False,
          console=__console__ , icon=imgpath+'pulsar.ico')

coll = COLLECT( exe,
               a.binaries,
               [('splash.bmp', imgpath+'splash.bmp', 'DATA'),
		('pulsar.png', imgpath+'pulsar.png', 'DATA')],
               strip=False,
               upx=False,
               name='pypulsar-'+__version__)
