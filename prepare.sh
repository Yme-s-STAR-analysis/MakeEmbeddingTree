#!/bin/bash
if [ ! $# -eq 1 ];then
    echo "Need target path."
else
    target=$1
    echo "target path is $1."
    rm Csubmit.xml
    sed -e "s|TARGET|$1|" .tpl.xml > Csubmit.xml

    mkdir -p $1/report
    mkdir -p $1/csh
    mkdir -p $1/err
    mkdir -p $1/list
    mkdir -p $1/log
    mkdir -p $1/out

    rm -r  $1/report/* 
    rm -r  $1/csh/* 
    rm -r  $1/err/* 
    rm -r  $1/list/* 
    rm -r  $1/log/* 
    rm -r  $1/out/*
    
    rm -r  libZip*.zip
    rm -r  libZip*.package
    rm -r  *dataset
    rm -r  *session.xml
fi
