#!MC 1410
$!VarSet |MFBD| = '/home/pritam/Desktop/Github/CUDA-aware-MPI/Output'
$!VarSet |T| = 0
$!LOOP 161
$!PICK SETMOUSEMODE
  MOUSEMODE = SELECT
$!PAGE NAME = 'Untitled'
$!PAGECONTROL CREATE
$!PICK SETMOUSEMODE
  MOUSEMODE = SELECT
$!OPENLAYOUT  "|MFBD|/TemplateDendriticGrowth.lpk"
$!READDATASET  '"|MFBD|/Field-|T|.tec" '
  READDATAOPTION = NEW
  RESETSTYLE = NO
  VARLOADMODE = BYNAME
  ASSIGNSTRANDIDS = YES
  VARNAMELIST = '"X" "Y" "phi" "T"'
$!PRINTSETUP PALETTE = COLOR
$!EXPORTSETUP IMAGEWIDTH = 699
$!IF |T| <= 9
  $!EXPORTSETUP EXPORTFNAME = '|MFBD|/Frame-000|T|.png'
$!ELSE
  $!IF |T| <= 99
    $!EXPORTSETUP EXPORTFNAME = '|MFBD|/Frame-00|T|.png'
  $!ELSE
    $!EXPORTSETUP EXPORTFNAME = '|MFBD|/Frame-0|T|.png'
  $!ENDIF
$!ENDIF
$!EXPORT 
  EXPORTREGION = CURRENTFRAME
$!VarSet |T| += 1
$!ENDLOOP
$!RemoveVar |MFBD|