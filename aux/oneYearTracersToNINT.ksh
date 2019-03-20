#!/bin/ksh
# Please see usage notes below. For an example script to run this in a loop,
# see aux/runTracersToNINT_example.ksh
##################################################################################
# DUST HARDCODES:
# Break Clay tracer into 4 dust NINT inputs, by these fractions:
set -A clayCoeff --  0.009 0.081 0.234 0.676
# redust/rodust formerly read from existing DUST NINT file (was hardcoded anyway):
nsized=7
set -A redust -- 0.132 0.23 0.416 0.766 1.386 2.773 5.545
set -A rodust -- 2.5 2.5 2.5 2.5 2.65 2.65 2.65
##################################################################################

function usage {
echo " "
echo " Prerequisites: * Monthly acc files are in the input subdirectory (averaged over years if needed)."
echo '                * NCO(perators) must be installed in the ${NCO} location passed to this script.'
echo " "
echo " Pass this script arguments:"
echo "  1: the run name"
echo "  2: the string for averaged years (e.g. 1852-1858 or 2005)"
echo "  3: representativeYear (for example you might run a 'year 2000' run for"
echo "     model years 1852-1858), use 2000."
echo "  4: input directory (where you have the acc files and some new files will be"
echo "     created and destroyed in there as well)"
echo "  5: output directory where new NINT input files will be created"
echo "  6: a reference year to measure time from"
echo "  7: the accelaration due to gravity (for air mass calculations; made argument"
echo "     in case that helps with other planets.)"
echo "  8: the bin(aries) directory to use for the netCDF operators"
echo "  9: skip_ozone: > 0 means ozone creation will be skipped; < 0 means chop off that"
echo "     many layers from the top of the tracer ozone. For example to avoid conflict of"
echo "     layer edges with the rad code level input."
echo " 10: skip_aerosols: non-zero means aerosols will be skipped"
echo " 11: skip_BCalbedo: non-zero means BC albedo calculation will be skipped"
echo " "
echo " Example:"
echo " $1 E4TcadiF40 1852-1858 2000 DATA NEW 1750 9.80665E0 /usr/local/other/SLES11.1/nco/4.4.4/intel-12.1.0.233/bin 0 1 1"
echo " will use files like DATA/SEP1852-1858.accE4TcadiF40.nc which represent a run"
echo " with year 2000 conditions, measuring time from Jan 1750, and make NINT input"
echo " files like:"
echo "    NEW/SUL_E4TcadiF40_kg_m2_144x90x40/2000.nc"
echo " and other aerosols (in this example, skipping ozone and BCalbedo). Using netCDF operators in:"
echo " /usr/local/other/SLES11.1/nco/4.4.4/intel-12.1.0.233/bin and acceleration due to"
echo " gravity of 9.80665E0 m s-2."
echo ""
echo " IMPORTANT NOTE: The ozone file that is created at most up to the top of the chemistry and hence is"
echo " intended for use only as an O3file2 with an O3file still listed with ozone information above this top."
echo " Furthermore, this functionality is currently only available on tracersE3 branch and downstream of that."
echo " To remind user of this, the created O3 filename gets prepended with 'E3_only' for now."
echo ""
echo " See runTracersToNINT_example.ksh as an example of how to run this script in a loop over years."
echo ""
echo " Feel free to contact Gregory.S.Faluvegi@nasa.gov."
echo " "
}

#---------------------------------------- functions ------------------------------------
function exitWithMessage {
  echo "Exiting because: $1" ; exit 1
}
function check {  # exits if there is more than one matching line for the grep:
  if [[ $( ncdump -h $1 | grep "$2" | cut -d= -f 2 | cut -d';' -f 1 | wc -l ) -ne 1 ]] ; then
    exitWithMessage 'Ambiguous grep'
  fi
}
function scrape { # gathers data from between the = sign and the semicolon for a grepped string:
  check $1 "$2"
  echo $( ncdump -h $1 | grep "$2" | cut -d= -f 2 | cut -d';' -f 1 )
}
function store { # stores the output of function scrape into a netCDF file:
  $NCO/ncap2 -A -s "${1}={${2}} ;" $3 $3
}
function removeIfExists {
  if [[ -f $1 ]] ; then rm $1; fi
}
function getPowers { # determines the power of 10 in units string - at least for aerosol mass:
  echo "1.E$( ncdump -h $1 | grep "${2}:units" | cut -d^ -f 2 | cut -d' ' -f 1 )"
}
function nocommas { # turns csv list into space-delimited:
  print -- $( echo $1 | sed 's/,/ /g' )
}

function buildAerosolNcapString {
  # Basically builds up the linear combinations of tracers needed from Susann'e program
  # including units conversion. Note it uses the main script variables!
  ncap2_string=''
  c=0
  for v in $( nocommas $needVars ) ; do
    tens=$( echo $( getPowers $taijl $v ) )
    if [[ $c -eq 0 ]] ; then
      ncap2_string="accum=${v}*${tens};"
    else
      ncap2_string="${ncap2_string} accum=accum+${v}*${tens};"
    fi
    let c+=1
  done
  # then convert units from (kg tracer / kg air) to (kg tracer / m2):
  ncap2_string="${ncap2_string} ${spec}=float(accum*airmass);"
}
function addDegenerateDimension {
  # Add a degenerate dimension to the variables. This is based on:
  # https://sourceforge.net/p/nco/discussion/9830/thread/0851525d/
  # function arg 1 = file, arg 2 = variable, arg 3 = temporary file
  mv $1 $3 ; $NCO/ncwa -a $2 $3 $1 ; $NCO/ncecat -O -u $2 $1 $1 ; rm $3
}

#--------------------- parse user args and check/prep in/out dirs -----------------------
echo "setup...\c"
if [[ $# -ne 11 ]] ; then usage $0 ; exitWithMessage "Expecting 11 arguments." ; fi
run=$1 ; avgYears=$2 ; repYear=$3 ; userInDir=$4 ; outDir=$5 
refYear=$6 ; grav=$7 ; NCO=$8 ; skipO3=$9 ; skipAero=${10} ; skipAlb=${11}

targets=''
if [[ $skipO3 -le 0 ]] ; then
  targets="$targets O3"
fi
if [[ $skipAero -eq 0 ]] ; then
  targets="$targets SSA SUL dust0 dust1 dust2 dust3 dust4 dust5 dust6 NIT OCA BCB BCA"
fi
if [[ $skipAlb  -eq 0 ]] ; then
  targets="$targets BCdalbsn"
fi

# Including dealing with whether user is using absolute or relative paths:
if [[ $( dirname $outDir ) != '.' ]] ; then
  echo "The script would rather create output locally. Is it OK to put output in:"
  echo "./$( basename $outDir ), instead of $outDir ?"
  echo "[Press return to continue or control-C to quit]" ; read nothing
fi
outDir=$( basename $outDir ) ; if [[ ! -d ${outDir} ]] ; then mkdir ${outDir} ; fi

base=${PWD}
inDir=___working
if [[ -d ${inDir} ]] ; then
  exitWithMessage "The directory ./${inDir} exists. Please remove it and try again."
else
  mkdir $inDir
  if [[ -d ./$userInDir ]] ; then  # data dir is local
    cd $inDir; ln -s ${base}/${userInDir}/???${avgYears}.{acc,ajl,taijl,aij,taij}${run}.nc . ; cd $base
  else
    cd $inDir; ln -s ${userInDir}/???${avgYears}.{acc,ajl,taijl,aij,taij}${run}.nc . ; cd $base
  fi
fi

echo "=-=-=-=- Doing: ${run} ${avgYears} ('$repYear') -=-=-=-="
# Check whether to use the himem or normal scaleacc in case scaled files are needed:
scaler=scaleacc ; if [[ $( hostname ) == borg* ]] ; then scaler=${scaler}_himem ; fi

#--------------------- Gather vertical grid info from JAN acc --------------------------
echo "extract grid info...\c"
fz=_zgrid.nc ; removeIfExists $fz
acc=${inDir}/JAN${avgYears}.acc${run}.nc
# Determine number of layers, then create a netCDF file from scratch:
plm=$( scrape $acc 'iparam:lm =' ) ; plmm1=$( echo "$plm - 1" | bc ) ; plmp1=$( echo "$plm + 1" | bc )
ls1=$( scrape $acc 'iparam:ls1 =' ) ; ls1m1=$( echo "$ls1 - 1" | bc )
echo "netcdf $( basename $fz .nc){\n dimensions: plm = ${plm};\n ple = ${plmp1};\n variables:\n double mfix(plm);\n double mfrac(plm);\n double plbot(ple);\n}" > $( basename $fz .nc).cdl
ncgen -b -o $fz $( basename $fz .nc).cdl ; rm $( basename $fz .nc).cdl
# Read the layering info from the acc file real parameter database, and store in netCDF for manipulation:
for v in mfix mfrac plbot ; do store $v "$( scrape $acc $v)" $fz ; done
$NCO/ncap2 -A -s "mtop=plbot($plm)*(100./$grav);" $fz $fz
$NCO/ncap2 -A -s "mfixs=plbot($ls1m1)*(100./$grav)-mtop;" $fz $fz
$NCO/ncatted -a units,plbot,o,c,'millibar' $fz
# Determine from whether or not first mfix element is negative whether this is a standard hybrid case:
stdhyb=1 ; [[ $( echo $( scrape $acc mfix ) | cut -d, -f 1 ) -lt 0 ]] || stdhyb=0

#--------------------- Main tracer work ------------------------------------------------
echo "tracers work...."
fn=_new.nc ; removeIfExists $fn
# Loop over months, concatting tracers into a 12-month file:
m=0
for mon in JAN FEB MAR APR MAY JUN JUL AUG SEP OCT NOV DEC ; do
  echo "${mon}: \c"
  fm=__model_${mon}
  for f in ${fm}.nc ${fm}.cdl __temp.nc ; do removeIfExists $f ; done
  # Make a new netcdf file for current month:
  echo "netcdf $fm {\ndimensions: time = UNLIMITED;\nvariables:\nfloat time(time);\ndata:\ntime=0.;\n}" > ${fm}.cdl
  ncgen -b -o ${fm}.nc ${fm}.cdl ; rm ${fm}.cdl
  # Give time the correct units:
  $NCO/ncatted -a units,time,o,c,"months since ${refYear}-01" ${fm}.nc

  # Create the scaled files, WARNING: ONLY CREATES IF THEY DON'T EXIST ALREADY!
  acc=${inDir}/${mon}${avgYears}.acc${run}.nc
  taijl=${inDir}/${mon}${avgYears}.taijl${run}.nc
  taij=${inDir}/${mon}${avgYears}.taij${run}.nc
  aij=${inDir}/${mon}${avgYears}.aij${run}.nc
  ajl=${inDir}/${mon}${avgYears}.ajl${run}.nc
  if [[ ! -f $taijl ]] ; then $scaler $acc taijl ; mv $( basename $taijl ) $taijl ; fi
  if [[ ! -f $aij ]] ; then $scaler $acc aij ; mv $( basename $aij ) $aij ; fi
  if [[ ! -f $ajl ]] ; then $scaler $acc ajl ; mv $( basename $ajl ) $ajl ; fi
  if [[ $skipAlb -eq 0 ]] ; then
    if [[ ! -f $taij ]] ; then $scaler $acc taij ; mv $( basename $taij ) $taij ; fi
  fi
  if [[ $mon == 'JAN' ]] ; then
    $NCO/ncks -A -v level $taijl $fz
    $NCO/ncks -A -v plm $ajl $fz
    nlev=$( echo $( scrape $taijl 'level =' ) | tr -d ' ' )
    nlon=$( echo $( scrape $taijl 'lon =' ) | tr -d ' ' )
    nlat=$( echo $( scrape $taijl 'lat =' ) | tr -d ' ' )
  fi
  if [[ $skipAero -eq 0 ]] ; then
    # Calculate or extract this month's 3D air mass (in kg m-2 per layer) for aerosol unit conversions:
    airmassExists=$( ncdump -h $taijl | grep 'float airmass' | wc -l )
    if [[ $airmassExists -ne 1 ]] ; then airmassExists=0 ; fi
    if [[ $airmassExists -eq 0 ]] ; then
      echo "airmass(make)...\c"
      # See explanations in model AtmLayering.F90 and ATMDYN_COM.F90.
      $NCO/ncks -A -v prsurf $aij $fz
      $NCO/ncap2 -A -s "airmass[\$level,\$lat,\$lon]=0.;" $fz $fz
      L=0
      while [[ $L -lt $plm ]]; do # size of level dim must be same as size of plm dim...
        if [[ $stdhyb -eq 1 ]] ; then
          $NCO/ncap2 -A -s "airmass($L,:,:)=mfix($L)+mfrac($L)*prsurf(:,:)*(100./$grav) ;" $fz $fz
        else
          $NCO/ncap2 -A -s "airmass($L,:,:)=mfix($L)+mfrac($L)*(prsurf(:,:)*(100./$grav)-mfixs-mtop) ;" $fz $fz
        fi
        let L+=1
      done
    else
      echo "airmass(get)...\c"
      $NCO/ncks -A -v airmass $taijl $fz
    fi
    # Put air mass into current monthly file:
    $NCO/ncks -A -v airmass $fz ${fm}.nc
  fi

  for spec in $targets ; do # loop over target species
    echo "${spec}...\c"
    units='kg/m2'
    case $spec in
    O3)
      needVars='O3_cm_atm' ; long_name='ozone amount' ; units='cm-atm'
      ncap2_string="${spec}=O3_cm_atm ;" ;;
    SUL)
      needVars='SO4' ; long_name='Sulfate mass' ;;
    NIT)
      needVars='NO3p' ;  long_name='Nitrate mass' ;;
    BCB) # careful as this one alters the extracted model var b/c same name as target.
      needVars='BCB' ; long_name='Black Carbon Bio Mass mass' ;;
    BCA)
      needVars='BCII,BCIA' ; long_name='Black Carbon Fossil Fuel mass' ;;
    OCA)
      needVars='OCII,OCIA,OCB,apinp1a,apinp2a,isopp1a,isopp2a' ; long_name='Organic Aerosol mass' ;;
    SSA)
      needVars='seasalt1,seasalt2' ; long_name='Sea Salt mass' ;;
    dust[0-3])
      needVars='Clay'
      tens=$( echo $( getPowers $taijl Clay ) )
      ncap2_string="${spec}=float(${clayCoeff[${spec:4:1}]}*${tens}*Clay*airmass);" ;;
    dust[4-6])
      needVars="Silt$( echo "${spec:4:1} - 3" | bc )" ;;
    BCdalbsn)
      needVars="alb_BC,sunlit_snow_freq" ;;
    esac
    case $spec in
      dust*) ; long_name='Dust mass' ;;
    esac
    case $spec in
      SUL|NIT|BCB|BCA|OCA|SSA|dust4|dust5|dust6) ; buildAerosolNcapString ;;
    esac

    # Gather needed variables and create the rad code input variable:
    case $spec in
    BCdalbsn)
      $NCO/ncks -A -v $needVars $taij ${fm}.nc
      $NCO/ncatted -a missing_value,alb_BC,d,, ${fm}.nc
      $NCO/ncrename -v alb_BC,${spec}_orig ${fm}.nc
      $NCO/ncap2 -A -s "${spec}=${spec}_orig;" ${fm}.nc ${fm}.nc
      ;;
    *)
      $NCO/ncks -A -v $needVars $taijl ${fm}.nc
      $NCO/ncap2 -A -s "$ncap2_string" ${fm}.nc ${fm}.nc
      $NCO/ncatted -a long_name,${spec},o,c,"$long_name" ${fm}.nc
      $NCO/ncatted -a units,${spec},o,c,"$units" ${fm}.nc
      ;;
    esac
  done

  if [[ $skipAero -eq 0 ]] ; then
    # remove airmass and surface pressure from files so they can be redefined next month:
    mv $fz __temp.nc ; $NCO/ncks -x -v prsurf,airmass __temp.nc $fz ; rm __temp.nc
    mv ${fm}.nc __temp.nc ; $NCO/ncks -x -v airmass __temp.nc ${fm}.nc ; rm __temp.nc
  fi

  echo "concat..."
  addDegenerateDimension ${fm}.nc time __temp.nc
  # Give time a value measured from reference year's Jan:
  m3=$( echo " 12 * ( $repYear - $refYear ) + $m " | bc )
  $NCO/ncap2 -A -s "time=${m3} ;" ${fm}.nc ${fm}.nc
  # Build up a 12-month file; removing current month's file:
  if [[ -f $fn ]] ; then  # after first time, append to the file...
    $NCO/ncrcat $fn ${fm}.nc __temp.nc ; mv -f __temp.nc $fn
  else                    # first time, create it...
    $NCO/ncks ${fm}.nc $fn
  fi
  rm ${fm}.nc
  let m+=1
done # months loop


#------------------ Extract species-specific NINT input --------------------------------
# Section could be greatly simplified if we harmonize formats... I wanted to match existing files:
echo "output..."
for spec in $targets ; do # loop over target species
  case $spec in
  O3)
    topMidIndex=${nlev} # e.g. starting at 1
    if [[ $skipO3 -lt 0 ]] ; then
      topMidIndex=$( echo " ${nlev} ${skipO3} " | bc )
      topEdgeIndex=$( echo " ${topMidIndex} + 1 " | bc )
    fi
    ds=${outDir}/E3_only_${spec}_${run}_cm-atm_${nlon}x${nlat}x${topMidIndex}
    ;;
  SUL|NIT|OCA|BCB|BCA|SSA|dust*) ; ds=${outDir}/${spec}_${run}_kg_m2_${nlon}x${nlat}x${nlev} ;;
  BCdalbsn) ; ds=${outDir}/${spec}_${run}_percent_${nlon}x${nlat} ;;
  esac
  fs=${ds}_${repYear}.nc
  $NCO/ncks -v $spec $fn $fs
  $NCO/ncatted -a cell_methods,time,d,, $fs
  case $spec in
  O3) # =========== special to O3 file ==============
    $NCO/ncatted -a long_name,lat,o,c,'lat' $fs
    $NCO/ncatted -a long_name,lon,o,c,'lon' $fs
    $NCO/ncrename -d level,plm $fs
    # Add plbot from the vertical info file and change its name to ple because it will be
    # read as that from the rad code and possibly used for interpolation:
    $NCO/ncks -A -v plbot $fz $fs
    $NCO/ncrename -v plbot,ple $fs
    $NCO/ncatted -a long_name,ple,o,c,'layer edge pressure' $fs
    # remove level variable:
    mv $fs __temp.nc ; $NCO/ncks -x -v level __temp.nc $fs ; rm __temp.nc
    # Option to chop off some layers of ozone tracer from top of file:
    if [[ $skipO3 -lt 0 ]] ; then
      echo "but only use first ${topMidIndex} levels for ozone..." 
      mv $fs __temp.nc
      # next line use (-F) fortran indexing:
      $NCO/ncks -O -F -d plm,1,${topMidIndex} -d ple,1,${topEdgeIndex} __temp.nc $fs  
      rm __temp.nc
    fi
    ;;
  SUL|NIT|OCA|BCB|BCA|SSA|dust*) # ===== special to aerosol files =====
    $NCO/ncatted -a long_name,time,o,c,"months since ${refYear}-01" $fs
    $NCO/ncatted -a long_name,lat,o,c,'Latitude' $fs
    $NCO/ncatted -a long_name,lon,o,c,'Longitude' $fs
    # rename level dimension and variable to lev and fill with plm values; remove plm.
    $NCO/ncrename -d level,lev $fs
    $NCO/ncrename -v level,lev $fs
    $NCO/ncks -A -v plm $fz $fs
    $NCO/ncap2 -A -s "lev=plm;" $fs $fs
    $NCO/ncatted -a units,lev,o,c,'millibar' $fs
    $NCO/ncatted -a long_name,lev,o,c,'Level' $fs
    mv $fs __temp.nc ; $NCO/ncks -x -v plm __temp.nc $fs ; rm __temp.nc
    case $spec in
    SSA) # ===== special to sea salt file =====
      # Add plbot as plbaer. Formerly read from exising SSA NINT file.
      $NCO/ncks -A -v plbot $fz $fs
      $NCO/ncap2 -A -s 'plbaer=float(plbot);' $fs $fs
      $NCO/ncatted -a long_name,plbaer,o,c,'pressure level' $fs
      $NCO/ncrename -d ple,levp1 $fs
      mv $fs __temp.nc ; $NCO/ncks -x -v plbot __temp.nc $fs ; rm __temp.nc
      ;;
    dust*) # ===== special to dust file =====
      # Add plbot as plbdust. Formerly read from exising DUST NINT file.
      $NCO/ncks -A -v plbot $fz $fs
      $NCO/ncap2 -A -s 'plbdust=float(plbot);' $fs $fs
      $NCO/ncatted -a long_name,plbdust,o,c,'layer edge pressure' $fs
      $NCO/ncrename -d ple,plbdust $fs
      # fix the time dimension, create size dimension and size variables, filling those for the current dust bin:
      $NCO/ncks -O --fix_rec_dmn time $fs $fs
      $NCO/ncap2 -A -s 'defdim("nsized",1);' $fs $fs
      # move below addDegenerateDimension $fs nsized __temp.nc
      $NCO/ncap2 -A -s "redust[\$nsized]=float(0.); rodust[\$nsized]=float(0.); nsized[\$nsized]=float(0.)" $fs $fs
      i0=${spec:4:1} ; i1=$( echo "$i0 + 1" | bc)
      $NCO/ncap2 -A -s "redust=${redust[${i0}]} ;" $fs $fs
      $NCO/ncap2 -A -s "rodust=${rodust[${i0}]} ;" $fs $fs
      $NCO/ncap2 -A -s "nsized=${i1} ;" $fs $fs
      $NCO/ncatted -a units,redust,o,c,' ' $fs ; $NCO/ncatted -a units,rodust,o,c,' ' $fs
      $NCO/ncatted -a long_name,redust,o,c,redust $fs ; $NCO/ncatted -a long_name,rodust,o,c,rodust $fs
      addDegenerateDimension $fs nsized __temp.nc
      # make dust size the outer dimension:
      $NCO/ncpdq -A -v ${spec} -a nsized,time,lev,lat,lon $fs $fs
      mv $fs __temp.nc ; $NCO/ncks -x -v plbot __temp.nc $fs ; rm __temp.nc
      ;;
    esac
    ;;
  BCdalbsn) # ===== special to BC albedo =====
    $NCO/ncks -A -v ${spec}_orig,sunlit_snow_freq $fn $fs
    # obtain the location of the TOPO file from the JAN acc file, extract glacial ice fraction:
    acc=${inDir}/JAN${avgYears}.acc${run}.nc
    topo=$( echo $( scrape $acc '_file_topo' ) | cut -d\" -f 2 )
    $NCO/ncks -A -v fgice $topo $fs
    # Repair the lat variable (tradiation has 90 and -90 as the first and last values;
    # importing those topo variables in last line changes those to 89 and -89):
    $NCO/ncks -A -v lat $taij $fs
    echo "Filling missing $spec data. Very slow as programmed. Expect about 3 minutes:"
    # Write a ncap script mimicking Max Kelley's fillit FORTRAN routine. Basically,
    # we are filling missing BCdalbsn points by linear interpolation from nearest bracketing
    # non-missing points in the N-S direction (or from 0 if no more non-missing available on
    # one side). But then for points that needed filling but have glacial ice present, reset to 0:
    maxmiss=-1.e20 ; externalScript=_fillit.nco
    echo "for(*k=0;k<12;k++){"                                                               > $externalScript
    echo " for(*i=0;i<${nlon};i++){"                                                        >> $externalScript
    echo "  while(${spec}(k,:,i).min()<${maxmiss}){"                                        >> $externalScript
    echo "    for(*j1=0;j1<${nlat}-1;j1++){"                                                >> $externalScript
    echo "      if(${spec}(k,j1,i)<${maxmiss}) break;"                                      >> $externalScript
    echo "    }"                                                                            >> $externalScript
    echo "    if(j1==0){*v1=0.0;} else{j1-=1; *v1=${spec}(k,j1,i);}"                        >> $externalScript
    echo "    for(*j2=j1+1;j2<${nlat};j2++){"                                               >> $externalScript
    echo "      if(${spec}(k,j2,i)>=${maxmiss}||j2==${nlat}-1) break;"                      >> $externalScript
    echo "    }"                                                                            >> $externalScript
    echo "    if(${spec}(k,j2,i)<${maxmiss}){*v2=0.0;} else{*v2=${spec}(k,j2,i);}"          >> $externalScript
    echo "    for(*j=j1;j<j2+1;j++){"                                                       >> $externalScript
    echo "      *wt1=(j2-j).float()/(j2-j1).float();"                                       >> $externalScript
    echo "      ${spec}(k,j,i)=wt1*v1+(1.-wt1)*v2;"                                         >> $externalScript
    echo "    }"                                                                            >> $externalScript
    echo "  }"                                                                              >> $externalScript
    echo "  *var_tmp=${spec}(k,:,i) ; *a=${spec}_orig(k,:,i) ; *b=fgice(:,i);"              >> $externalScript
    echo "  where (a < ${maxmiss} && b > 0.) var_tmp=0.;"                                   >> $externalScript
    echo "  ${spec}(k,:,i)=var_tmp;"                                                        >> $externalScript
    echo "  ram_delete(var_tmp,a,b);"                                                       >> $externalScript
    echo " }"                                                                               >> $externalScript
    echo "}"                                                                                >> $externalScript
    date
    $NCO/ncap2 -h -A -S $externalScript $fs $fs
    date
    rm $externalScript
    mv $fs __temp.nc ; $NCO/ncks -x -v fgice,${spec}_orig __temp.nc $fs ; rm __temp.nc
    ;;
  esac
  # Make the DIR/YYYY.nc timestream version:
  if [[ ! -d ${ds} ]] ; then mkdir ${ds} ; fi
  # _could_ change mv to cp on next line to retain (single year) NON-dir/YYYY.nc version:
  [[ ! -e ${ds}/${repYear}.nc ]] || exitWithMessage "${ds}/${repYear}.nc already exists."
  mv $fs ${ds}/${repYear}.nc
  if [[ ${spec:0:4} != 'dust' ]] ; then
    ls -l ${ds}/${repYear}.nc  # show user the file
  else
    echo "${spec}..."
  fi

done # targets

#------------------ Combine dusts into one file ----------------------------------------
if [[ $skipAero -eq 0 ]] ; then
  echo "\n combining dusts..."
  dd=${outDir}/DUST_${run}_kg_m2_${nlon}x${nlat}x${nlev}
  if [[ ! -d $dd ]] ; then mkdir $dd ; fi
  fd=${dd}/${repYear}.nc
  flist=''
  for dustX in dust0 dust1 dust2 dust3 dust4 dust5 dust6 ;do
    fx=${outDir}/${dustX}_${run}_kg_m2_${nlon}x${nlat}x${nlev}/${repYear}.nc
    flist="$flist $fx"
    $NCO/ncrename -v ${dustX},DUST $fx # common name for ncrcat below.
  done
  $NCO/ncrcat $flist $fd
  # return time to be the record dim in final file and correct dimensions order:
  $NCO/ncks -O --fix_rec_dmn nsized $fd $fd
  # next line fails if in.nc and out.nc are the same file, so separate DUST var...
  $NCO/ncpdq -A -v DUST -a time,nsized,lev,lat,lon $fd __temp.nc
  # ... and non-dust vars:
  $NCO/ncks -x -v DUST $fd __rest.nc
  rm $fd
  # ... and recombine them:
  $NCO/ncks -A __rest.nc $fd
  $NCO/ncks -A __temp.nc $fd
  $NCO/ncks -O --mk_rec time $fd $fd
  $NCO/ncrename -v nsized,bins $fd
  $NCO/ncatted -a units,bins,o,c,'bins' $fd
  $NCO/ncatted -a long_name,bins,o,c,'dust size bins' $fd
  $NCO/ncatted -a cell_methods,redust,d,, $fd ; $NCO/ncatted -a cell_methods,rodust,d,, $fd
  $NCO/ncatted -a cell_methods,bins,d,, $fd
  rm __temp.nc __rest.nc
  ls -l $fd
  for dustX in dust0 dust1 dust2 dust3 dust4 dust5 dust6 ;do
    fx=${outDir}/${dustX}_${run}_kg_m2_${nlon}x${nlat}x${nlev}/${repYear}.nc
    rm $fx ; rmdir ${outDir}/${dustX}_${run}_kg_m2_${nlon}x${nlat}x${nlev}
  done
fi

#--------------------------- Clean-up --------------------------------------------------
echo "cleanup..."
rm $fz $fn
cd $inDir ; rm -f ???${avgYears}.{aij,ajl,taijl,taij,acc}${run}.nc ; cd $base
rmdir ${inDir}
