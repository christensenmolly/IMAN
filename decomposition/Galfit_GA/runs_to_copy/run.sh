#!/bin/bash

gnome-terminal --tab -x bash -c "sh run_1.sh $1" &
gnome-terminal --tab -x bash -c "sh run_2.sh $1" &
gnome-terminal --tab -x bash -c "sh run_3.sh $1" &
gnome-terminal --tab -x bash -c "sh run_4.sh $1" &

