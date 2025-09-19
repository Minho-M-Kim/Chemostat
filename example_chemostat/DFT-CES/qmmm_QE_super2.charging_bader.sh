#!/bin/bash

# Global variables
DO_PARALLEL="mpirun"
QMTYPE="scf" # scf or opt
QMMAXSTEP="5" # QM max trial for "opt"
QMMMMAXSTEP="4" # no of QMMM loop
SUPERCELL=(11 11 1) # supercell factor x, y, z
DIPOLEDIR="3" # dipole correction direction x=1 y=2 z=3
DIPOLEPOS="0.80" # dipole correction fractional position along DIPOLEDIR
LAMMPS="(Your Own LAMMPS path)"
QEPW="(Your Own QE pw.x path)"
QEPP="(Your Own QE pp.x path)"
CUBESUB="(Your Own cube_sub path)"
CUBEADD="(Your Own cube_add path)"
CUBEMULTI="(Your Own cube_multi path)"
MDIPC="(Your Own cube_multi path)"
CHG2POT="(Your Own chg2pot path)"
LAMMPSIN="base.in.lammps"
QMIN="base.pw.in"
QMIN2="base.pp.in"
LAMMPSDATA="data.temp"
LAMMPSRESTART="Graphene_1000H2O_3cation.NVT.5000000.restart"
MDEQUIL="2000000" # for pzc, set to 500000 also average step 
MDAVERAGE="2000000" # they are multiplied in loop $maxavgstep times
maxavgstep="5" #avgstep after mm_0
maxavgstep2="5" #avgstep after mm_3
NSOLVATOM="18194" # Number of solvent atom (lammps id: 1..N) + ions
TOTCHG="0.0" # (-):put more electrons, for QM cal , should be "0.0" for PZC
CONCEN="dilute2" # should write dilute if only compensation ions are contained... 
mpclayer="11.75" #bohr unit mpc layer
totlayer="4" 
mpcone="1"
adsorbate="0"
initialqm=0 #initial qm exist yes=1 no=0
finalmm=0 #final md do or not yes=1 no=0
skipeq=0
startavg=1
qmmmstep=0
firstrun=0
SOFT="0.532999"

for ((qmmmstep=$qmmmstep; qmmmstep<$QMMMMAXSTEP; qmmmstep++))
do
  echo "######### Starting $qmmmstep QMMM step #########"$'\n'
  # Make saving directory
  mkdir qm_$qmmmstep
  mkdir mm_$qmmmstep

  if [ $firstrun -eq 1 ]; then
    # Parsing QE input info
    natoms_qm=`grep nat $QMIN | sed 's/,//' | awk '{print $3}'`
    echo "### Number of atom in QE input: $natoms_qm"$'\n'
    atoms_qm=(`awk -v naqm=$natoms_qm 'BEGIN{i=1}{if($1=="ATOMIC_POSITIONS") n=NR; if(n!=0&&NR>n&&NR<n+naqm+1){if($5==""){print $0,"1 1 1"}else{print $0}}}' $QMIN`) # QE input format
    solute=(`awk -v namm=$NSOLVATOM -v naqm=$natoms_qm 'BEGIN{i=1}{if($1=="ATOMIC_POSITIONS") n=NR; if(n!=0&&NR>n&&NR<n+naqm+1){print namm+i,$1,$2,$3,$4,0;i++}}' $QMIN`) # new_lammps_id, species, x, y, z
    qm_cells=(`awk '{if($1=="CELL_PARAMETERS") n=NR; if(n!=0&&NR>n&&NR<=n+3) print $(NR-n)}' $QMIN`)
    
    # lammps restart to data
    #$DO_PARALLEL $LAMMPS -r $LAMMPSRESTART data.temp

    awk -v namm=$NSOLVATOM '{if(NF==10&&$1>namm) print $0}' data.temp > data.solute
    awk -v namm=$NSOLVATOM '{if(!(NF==10&&$1>namm)) print $0}' data.temp > data.solvent

    # Parsing cell box
    cells=(`awk '{if(NR>=10&&NR<=12) printf "%.6f ", $2-$1}' data.temp`)
    echo "### Lammps cell parameters : ${cells[@]}"$'\n'

    # Finding lammps atom type for QM atoms
    for ((i=0; i<$natoms_qm; i++))
    do
      while read line
      do
        temp=""
        temp=`echo $line ${solute[((6*i+2))]} ${solute[((6*i+3))]} ${solute[((6*i+4))]} ${cells[@]} | awk '{dist=sqrt(($5+$8*$14-$11)**2+($6+$9*$15-$12)**2+($7+$10*$16-$13)**2);if(dist<0.1) print $3}'`
        if [ "$temp" != "" ]; then
          solute[((6*i+1))]=$temp
          break
        fi
      done < data.solute
    done

    # Finding lammps atom charge for QM atoms
    for ((i=0; i<$natoms_qm; i++))
    do
      while read line
      do
        temp=""
        temp=`echo $line ${solute[((6*i+2))]} ${solute[((6*i+3))]} ${solute[((6*i+4))]} ${cells[@]} | awk '{dist=sqrt(($5+$8*$14-$11)**2+($6+$9*$15-$12)**2+($7+$10*$16-$13)**2);if(dist<0.1) print $4}'`
        if [ "$temp" != "" ]; then
          solute[((6*i+5))]=$temp
          break
        fi
      done < data.solute
    done

    firstrun=0
  fi # 0th step initial parsing end

  # Make pw.in based on $QMIN 
  if [ "$QMTYPE" == "opt" ]; then
    awk -v ats="${atoms_qm[*]}" -v naqm=$natoms_qm 'BEGIN{split(ats,a);i=0}{if($1=="ATOMIC_POSITIONS") n=NR;if(n!=0&&NR>n&&NR<n+naqm+1){print a[7*i+1],a[7*i+2],a[7*i+3],a[7*i+4],a[7*i+5],a[7*i+6],a[7*i+7];i++}else{print $0}}' $QMIN > pw.in
  else
    cp $QMIN pw.in
  fi
  if [ $qmmmstep -gt 0 ]; then
    sed -i "s/.*\&CONTROL.*/&\ndft_ces = .true./" pw.in
    sed -i "s/.*\&CONTROL.*/&\nrho_ces = '.\/empty.cube'/" pw.in
    sed -i "s/.*\&CONTROL.*/&\nrho_ces_ion = '.\/MOBILE_final.cube'/" pw.in
    sed -i "s/.*\&CONTROL.*/&\nrho_ces_rep = '.\/empty.cube'/" pw.in
    sed -i "s/.*\&SYSTEM.*/&\ntot_charge = $TOTCHG/" pw.in
    #if [ $qmmmstep -gt 1 ]; then 
      #sed -i "s/.*\&ELECTRONS.*/&\nstartingpot = 'file'/" pw.in 
      #sed -i "s/startingwfc = 'atomic'/startingwfc = 'file'/" pw.in
      if [ "$QMTYPE" == "opt" ]; then
        sed -i "s/.*\&CONTROL.*/&\nforc_conv_thr = 1.0D-3/" pw.in
        sed -i "s/.*\&CONTROL.*/&\nnstep = 500/" pw.in
        sed -i "s/.*calculation.*/calculation = 'relax'/" pw.in
        sed -i "s/.*ATOMIC_SPECIES.*/\&IONS\n&/" pw.in
        sed -i "s/.*ATOMIC_SPECIES.*/\/\n&/" pw.in
      fi
    #fi
    cp ./mm_$((qmmmstep-1))/MOBILE_final.cube ./
    cp ./mm_$((qmmmstep-1))/empty.cube ./
   #cp ./mm_$((qmmmstep-1))/rep.cube ./
  fi
 
 # if [ "$TOTCHG" == "0.0" ]; then 
 #   if [ $qmmmstep -eq 1 ]; then 
 #     sed -i "s/.*\&ELECTRONS.*/&\nstartingpot = 'file'/" pw.in 
 #     sed -i "s/.*\&ELECTRONS.*/&\nstartingwfc = 'file'/" pw.in
 #   fi
 # fi

  # Run QM calculation
  qmcnt=1
  forcthr=1
  finished=""
  for ((i=0; i<$QMMAXSTEP; i++))
  do
    echo "### QM calculation starting No. $qmcnt/$QMMAXSTEP in $qmmmstep QMMM iterations"$'\n'
    if [ $initialqm -eq 0 ]; then
     $DO_PARALLEL $QEPW < pw.in > pw.out
     if [ "$QMTYPE" == "scf" ]; then
       finished=`grep "JOB DONE" pw.out`
       if [ "$finished" != "" ]; then 
         echo "### QM calculation done"$'\n'
         break
       else
         echo "### QM calculation aborted, check output file"$'\n'
         exit
       fi
     elif [ $qmmmstep -eq 0 ]; then
       finished=`grep "JOB DONE" pw.out`
       if [ "$finished" != "" ]; then 
         echo "### QM calculation done"$'\n'
         break
       else
         echo "### QM calculation aborted, check output file"$'\n'
         exit
       fi
     else
       finished=`grep "bfgs converged in" pw.out`
       if [ "$finished" == "" ]; then
         let forcthr=forcthr*2
         echo "### QM calculation has not been converged, increasing froc_conv_thr to ${forcthr}D-3"$'\n'
         sed -i "s/.*forc_conv_thr.*/forc_conv_thr = ${forcthr}D-3/" pw.in
       else
         echo "### QM calculation has been converged"$'\n'
         break
       fi
       if [ "$finished" == "" -a $qmcnt -eq $QMMAXSTEP ]; then
         echo "### QM calculation failed, aborting whole loop."$'\n'
         exit
       fi
       let qmcnt=qmcnt+1
     fi
    else
     echo "###QM calculation in $qmmmstep QMMM iterations already exists"
     break
    fi
  done
  if [ $initialqm -eq 0 ]; then
   cp v_saw.cube v_saw_pw.cube
   cp pw.in pw.out v_saw_pw.cube v_md.cube v_md_ion.cube v_hartree.cube qm_$qmmmstep
  fi
  # Generating QM solute potential and backup QM results
  echo "### Generating QM potential .. "$'\n'
  # makeing pp.in
    cp $QMIN2 pp.pot.in
    cp $QMIN2 pp.rho.in
    sed -i "s/.*plot_num.*/plot_num = 0/" pp.rho.in
    sed -i "s/.*fileout.*/fileout = 'val.rho.cube'/" pp.rho.in

  #
  if [ $initialqm -eq 0 ]; then
   if [ "$TOTCHG" == "0.0" ] || [ $qmmmstep -eq 0 ]; then
     $DO_PARALLEL $QEPP < pp.pot.in > pp.pot.out
     ppdone1=`ls solute.cube`
     $DO_PARALLEL $QEPP < pp.rho.in > pp.rho.out
     ppdone2=`ls val.rho.cube`
     if [ "$ppdone1" == "" ] || [ "$ppdone2" == "" ] ; then
       echo "### QM post processing failed, aborting whole qmmm loop"$'\n'
       exit
     else
       cp val.rho.cube refval.rho.cube 
       $CUBEADD solute.cube v_saw_pw.cube
       $CUBEADD add.cube v_md_ion.cube
       cp add.cube total.cube
       $CUBEADD solute.cube v_saw.cube
       cp add.cube solute.cube
       cp pp.pot.in pp.rho.in pp.pot.out pp.rho.out v_saw.cube solute.cube total.cube refval.rho.cube qm_$qmmmstep
     fi
   else
     $DO_PARALLEL $QEPP < pp.pot.in > pp.pot.out
     ppdone1=`ls solute.cube`
     $DO_PARALLEL $QEPP < pp.rho.in > pp.rho.out
     ppdone2=`ls val.rho.cube`
     if [ "$ppdone1" == "" ] || [ "$ppdone2" == "" ] ; then
       echo "### QM post processing failed, aborting whole qmmm loop"$'\n'
       exit
     else
       cp solute.cube v_b+h.cube
       $CUBESUB val.rho.cube qm_0/refval.rho.cube
       dipolegrid=`awk -v dir=$DIPOLEDIR -v pos=$DIPOLEPOS '{if(NR==(dir+3)) printf "%d", $1*pos}' solute.cube`
       echo "### dipole correction has been applied: dir=$DIPOLEDIR, grid=$dipolegrid"
       $MDIPC subtracted.cube $DIPOLEDIR $dipolegrid $mpclayer $totlayer $mpcone $adsorbate
       $CHG2POT mdipc.cube $DIPOLEDIR
       $CUBEMULTI pot.cube 2
       $CUBEADD multiplied.cube qm_0/solute.cube
       cp add.cube solute.cube
       $CUBEADD v_b+h.cube v_md_ion.cube
       $CUBEADD add.cube v_saw_pw.cube
       cp add.cube total.cube
       cp pp.pot.in pp.rho.in pp.pot.out pp.rho.out v_saw.cube solute.cube total.cube val.rho.cube v_b+h.cube qm_$qmmmstep
       rm -f val.rho.cube
     fi
   fi
  fi
  initialqm=0
  # Updating QM optimized geometry
   if [ $qmmmstep -gt 0 -a "$QMTYPE" == "opt" ]; then
    solutefinal=(`awk '{if($0=="Begin final coordinates") tag=1 ;if($0=="End final coordinates") tag=0; if(tag==1 && (NF==4||NF==7)) print $2,$3,$4 }' pw.out`)
    for ((i=0; i<$natoms_qm; i++))
    do
      atoms_qm[((7*i+1))]=${solutefinal[((3*i))]}
      atoms_qm[((7*i+2))]=${solutefinal[((3*i+1))]}
      atoms_qm[((7*i+3))]=${solutefinal[((3*i+2))]}
      solute[((6*i+2))]=${solutefinal[((3*i))]}
      solute[((6*i+3))]=${solutefinal[((3*i+1))]}
      solute[((6*i+4))]=${solutefinal[((3*i+2))]}
    done
    echo "### QM optimized geometry has been parsed "$'\n'
  fi

  ##########################################
  if [ $((qmmmstep+1)) -eq $QMMMMAXSTEP ] && [ $finalmm -eq 0 ]; then
  echo "### MD at $qmmmstep step is not calculated"
  exit
  fi
  ##########################################

  # Making LAMMPS INPUT : modifying restart, geometry
  cp $LAMMPSIN in.lammps
  if [ "$QMTYPE" == "scf" ]; then
    sed -i "s/.*read_data.*/read_data $LAMMPSDATA/" in.lammps
  else # for QMTYPE=opt
    rm -f data.solute
    cnt=0
    for ((i=0; i<${SUPERCELL[0]}; i++))
    do
      for ((j=0; j<${SUPERCELL[1]}; j++))
      do
        for ((k=0; k<${SUPERCELL[2]}; k++))
        do
          awk -v cnt=$cnt -v su="$i $j $k" -v ats="${solute[*]}" -v cells="${qm_cells[*]}" -v n=$natoms_qm 'BEGIN{split(su,s);split(ats,a);split(cells,c);for(i=0;i<n;i++){printf "%d 444 %d %15.8f %15.8f %15.8f %15.8f 0 0 0\n",a[6*i+1]+cnt*n,a[6*i+2],a[6*i+6],a[6*i+3]+c[1]*s[1],a[6*i+4]+c[2]*s[2],a[6*i+5]+c[3]*s[3]}}' >> data.solute
          let cnt+=1
        done
      done
    done
    data_solute="`cat data.solute`"
    nl=`grep -n Atoms data.solvent | awk -F: '{print $1+2}'`
    awk -v solute="$data_solute" -v nl=$nl '{if(NR==nl) print solute;print $0}' data.solvent > data.lammps

    rm -f data.solute_chgcorr
    net_chg=$(awk '{if(NF==10)sum+=$4}END{print sum}' data.lammps)
    atag=0

    while read id temp type charge x y z a b c 
    do
     if [ "$charge" != "0.00000000" ] && [ $atag -eq 0 ]; then 
      corr_chg=$(echo "$charge - $net_chg" | bc -l)
      echo $id $temp $type $corr_chg $x $y $z $a $b $c >> data.solute_chgcorr
      atag=1
     else 
      echo $id $temp $type $charge $x $y $z $a $b $c >> data.solute_chgcorr
     fi
    done < data.solute
    data_solute="`cat data.solute_chgcorr`"
    nl=`grep -n Atoms data.solvent | awk -F: '{print $1+2}'`
    awk -v solute="$data_solute" -v nl=$nl '{if(NR==nl) print solute;print $0}' data.solvent > data.lammps

    sed -i "s/.*read_restart.*/read_data data.lammps/" in.lammps
  fi    
  
  cp ./qm_$qmmmstep/solute.cube ./

  # RUN LAMMPS calculation
  if [ $skipeq -eq 0 ]; then
    echo "### Running LAMMPS(emxext) for equilibration $qmmmstep QMMM iterations"$'\n'
    sed -i "s/.*run.*/run\t\t$MDEQUIL/" in.lammps
    sed -i "s/SOFT_INI/$SOFT/g" in.lammps
    $DO_PARALLEL $LAMMPS -in in.lammps > lammps.equil.out 
    cp in.lammps in.lammps.equil
    LAMMPSRESTART=`ls -lrt *.restart | tail -n 1 | awk '{print $9}'`
    $DO_PARALLEL $LAMMPS -r $LAMMPSRESTART $LAMMPSDATA #by jay
    SOFT=`grep Chemostat lammps.equil.out | tail -n 1 | awk '{print $7}'`
  fi

  skipeq=0

  cp $LAMMPSIN in.lammps
  if [ "$QMTYPE" == "scf" ]; then
    sed -i "s/.*read_data.*/read_data $LAMMPSDATA/" in.lammps
  else
    #sed -i "s/.*read_data.*/read_restart $LAMMPSRESTART/" in.lammps
    sed -i "s/.*read_data.*/read_data $LAMMPSDATA/" in.lammps #by jay
  fi

  echo "### Running LAMMPS(emdext) for averaging solvent charge density $qmmmstep QMMM iterations"$'\n'
  sed -i "s/.*run.*/run\t\t$MDAVERAGE/" in.lammps
  sed -i "s/SOFT_INI/$SOFT/g" in.lammps  
  if [ $qmmmstep -eq 0 ]; then 
      $DO_PARALLEL $LAMMPS -in in.lammps > lammps.average.out
      SOFT=`grep Chemostat lammps.average.out | tail -n 1 | awk '{print $7}'`
      cp SOLVENT.cube SOLVENT_final.cube
      cp cat.cube cat_final.cube
      cp ani.cube ani_final.cube
      cp catfrac.cube catfrac_final.cube
      cp anifrac.cube anifrac_final.cube
  else 
    if [ "$TOTCHG" == "0.0" ] && [ "$CONCEN" == "dilute" ]; then 
      $DO_PARALLEL $LAMMPS -in in.lammps > lammps.average.out
      SOFT=`grep Chemostat lammps.average.out | tail -n 1 | awk '{print $7}'`
      cp SOLVENT.cube SOLVENT_final.cube
      cp cat.cube cat_final.cube
      cp ani.cube ani_final.cube
      cp catfrac.cube catfrac_final.cube
      cp anifrac.cube anifrac_final.cube
    else
      if [ $qmmmstep -gt 3 ]; then
        maxavgstep=$maxavgstep2;
      fi
      for ((avgstep=$startavg; avgstep<=$maxavgstep; avgstep++)) 
      do
        cp $LAMMPSIN in.lammps
        cp in.lammps in.lammps.step$avgstep
        sed -i "s/.*run.*/run\t\t$MDAVERAGE/" in.lammps.step$avgstep
        sed -i "s/SOFT_INI/$SOFT/g" in.lammps.step$avgstep
        sed -i "s/.*read_data.*/read_data $LAMMPSDATA/" in.lammps.step$avgstep
        sed -i "s/{sname}.emd/{sname}.step$avgstep.emd/g" in.lammps.step$avgstep
        
        $DO_PARALLEL $LAMMPS -in in.lammps.step$avgstep > lammps.average.$avgstep.out
        SOFT=`grep Chemostat lammps.average.$avgstep.out | tail -n 1 | awk '{print $7}'`
        cp SOLVENT.cube SOLVENT_step$avgstep.cube
	LAMMPSTRJ=`ls -lrt *.lammpstrj | tail -n 1 | awk '{print $9}'`
	cp $LAMMPSTRJ step$avgstep.lammpstrj
	cp cat.cube cat_step$avgstep.cube
        cp ani.cube ani_step$avgstep.cube
        cp catfrac.cube catfrac_step$avgstep.cube
        cp anifrac.cube anifrac_step$avgstep.cube
        LAMMPSRESTART=`ls -lrt *.restart | tail -n 1 | awk '{print $9}'`
        $DO_PARALLEL $LAMMPS -r $LAMMPSRESTART $LAMMPSDATA
      done

      cp SOLVENT_step1.cube walksolv.cube
      cp cat_step1.cube walkcat.cube
      cp ani_step1.cube walkani.cube
      cp catfrac_step1.cube walkcatfrac.cube
      cp anifrac_step1.cube walkanifrac.cube
      cp walksolv.cube SOLVENT_1.TIME.cube
      cp walkcat.cube cat_1.TIME.cube
      cp walkani.cube ani_1.TIME.cube
      cp walkcatfrac.cube catfrac_1.TIME.cube
      cp walkanifrac.cube anifrac_1.TIME.cube

      for ((avgstep=2; avgstep<=$maxavgstep; avgstep++))
      do
        $CUBEADD walksolv.cube SOLVENT_step$avgstep.cube
        multiplier=$(echo "1/$avgstep" | bc -l)
        $CUBEMULTI add.cube $multiplier
        mv add.cube walksolv.cube
        mv multiplied.cube SOLVENT_$avgstep.TIME.cube
       
        $CUBEADD walkcat.cube cat_step$avgstep.cube
        $CUBEMULTI add.cube $multiplier
        mv add.cube walkcat.cube
        mv multiplied.cube cat_$avgstep.TIME.cube

        $CUBEADD walkani.cube ani_step$avgstep.cube
        $CUBEMULTI add.cube $multiplier
        mv add.cube walkani.cube
        mv multiplied.cube ani_$avgstep.TIME.cube

        $CUBEADD walkcatfrac.cube catfrac_step$avgstep.cube
        $CUBEMULTI add.cube $multiplier
        mv add.cube walkcatfrac.cube
        mv multiplied.cube catfrac_$avgstep.TIME.cube

        $CUBEADD walkanifrac.cube anifrac_step$avgstep.cube
        $CUBEMULTI add.cube $multiplier
        mv add.cube walkanifrac.cube
        mv multiplied.cube anifrac_$avgstep.TIME.cube
      done
      cp SOLVENT_$maxavgstep.TIME.cube SOLVENT_final.cube
      cp cat_$maxavgstep.TIME.cube cat_final.cube      
      cp ani_$maxavgstep.TIME.cube ani_final.cube
      cp catfrac_$maxavgstep.TIME.cube catfrac_final.cube
      cp anifrac_$maxavgstep.TIME.cube anifrac_final.cube
    fi
  fi
  cp in.lammps in.lammps.average
  LAMMPSRESTART=`ls -lrt *.restart | tail -n 1 | awk '{print $9}'`
  lammpsout=`ls -lrt lammps.*.out | tail -n 1 | awk '{print $9}'`
  finished=`grep "Please see the log.cite" $lammpsout`
  if [ "$finished" != "" ]; then 
    echo "### MM calculation done"$'\n'
  else
    echo "### MM calculation aborted, check output file"$'\n'
    exit
  fi

  startavg=1

  rm -f log.lammps
  $CUBEADD SOLVENT_final.cube cat_final.cube
  $CUBEADD add.cube ani_final.cube
  $CUBEADD add.cube catfrac_final.cube
  $CUBEADD add.cube anifrac_final.cube
  cp add.cube MOBILE_final.cube
  $CUBEMULTI MOBILE_final.cube 0
  cp multiplied.cube empty.cube
  cp data.temp in.lammps.* empty.cube *_final.cube *step*.cube *TIME*.cube *.lammpstrj $LAMMPSRESTART lammps.*.out mm_$qmmmstep
done # qmmm looP
