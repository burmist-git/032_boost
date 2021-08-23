#!/bin/bash

function nm_sh {

    liblist=$1
    
    for filename in $liblist ; do
	#nm -CD --defined-only $filename | grep -H --label= $filename 'test'
	#nm -CD --defined-only $filename | grep 'test::'
	#ff=$(nm -CD --defined-only $filename)
	echo " "
	echo "----> $filename"
	nm -CD --defined-only $filename | awk '{if($2=="T"){{print }}}' | grep -v _fini | grep -v _init | nl
	#nm -gDC $filename | grep _
	echo " "
    done
    
}

#nm -gDC *.so
#nm -S *.so
#nm_sh "*.so"

#nm_sh "/home/burmist/miniconda3/envs/boost_python/lib/*.so"
#nm_sh "/usr/lib/x86_64-linux-gnu/libboost_python*"

nm_sh hello_ext.so
nm_sh log_gaussian_heap_ext.so
nm_sh extending.so
nm_sh gaussian.so
