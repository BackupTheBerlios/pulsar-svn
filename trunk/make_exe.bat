@echo off

rem force linking and compiling
@del MODULES\*.*
@del OBJECTS\*.*
@del RELEASE\pulsar.exe  
@del errs.cf

AUTOMAKE fig=pulsar.fig 

call AMTEMP
call AMTEMP > errs.cf
DEL AMTEMP.BAT

@del *.dep
@del RELEASE\*.map

pause

