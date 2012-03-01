import sfle

Import( "*" )

execfile( "SConscript_maaslin.py" )

for filePCL in Glob( sfle.d( fileDirInput, "*" + sfle.c_strSufPCL ) ):
	Default( MaAsLin( filePCL ) )
