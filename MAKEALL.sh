#!/bin/sh
#////////////////////////////////////////////////////!
# * Script file to compile all modules and libraries
#////////////////////////////////////////////////////!

# src0 :  independent modules
# src1 :  depends only on public libraries
# src2 :  depends both on src0 and src1

cwd=$(pwd)

makes()
{
  echo '----' ${1} '----'
  cd ${1}
  make -f Makefile clean
  make -f Makefile
  make -f Makefile install
  make -f Makefile clean
  cd ${cwd}
}

if [ ${1} = "clean" ]; then

  rm -rf mod/*.mod
  rm -rf lib/*.a


elif [ ${1} = "fpicclean" ]; then

  rm -rf fpicmod/*.mod
  rm -rf fpiclib/*.a


elif [ ${1} = "install" ]; then

  rm -rf mod/*.mod
  rm -rf lib/*.a

  # src0
  makes src0

  # src1
  makes src1/linalg

  # src2
  makes src2/anaflat
  makes src2/anafull
  makes src2/nldd

fi


