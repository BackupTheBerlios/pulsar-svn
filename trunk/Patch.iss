; Script generated by the Inno Setup Script Wizard.
; SEE THE DOCUMENTATION FOR DETAILS ON CREATING INNO SETUP SCRIPT FILES!

#define MyAppName "PATCH 20070207 - pyPulsar"
#define MyAppVerName "pyPulsar-0.1.rc1"

[Setup]
AppName={#MyAppName}
AppVerName={#MyAppVerName}
DefaultDirName={pf}\{#MyAppVerName}
AllowNoIcons=true
InfoBeforeFile=D:\program_files\Pulsar\trunk\PATCH.txt
OutputBaseFilename=Patch-20070207-{#MyAppVerName}
Compression=lzma
SolidCompression=true
ShowLanguageDialog=auto
LanguageDetectionMethod=locale
WindowVisible=true
BackColor=clSilver
WizardImageFile=D:\program_files\Pulsar\trunk\Images\installer.bmp
WizardSmallImageFile=D:\program_files\Pulsar\trunk\Images\installer-small.bmp
WizardImageStretch=true
VersionInfoVersion=0.1.32
RestartIfNeededByRun=false
CreateAppDir=true
DirExistsWarning=no

[Languages]

[Tasks]

[Files]
; NOTE: Don't use "Flags: ignoreversion" on any shared system files
Source: f95pulsar.pyd; DestDir: {app}; Flags: overwritereadonly onlyifdestfileexists promptifolder

