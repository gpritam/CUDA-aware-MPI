#!/bin/bash

t=0;

flag=1;

while [ "$flag" = "1" ];
do
	if ls $(printf "Field-%d-*" $t) >/dev/null 2>&1;
	then
		cat $(printf "Field-%d-*" $t) > $(printf "Field-%d.tec" $t)
		rm $(printf "Field-%d-*" $t)
		t=$((t + 1))
	else
		flag=0
	fi
done
