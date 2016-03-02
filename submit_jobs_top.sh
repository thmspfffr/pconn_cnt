#!/bin/sh

# Submits n jobs to the torque queing system

for i in {1..14}
do
  echo 'Start Job' $i
  qsub submit_jobs.sh
  sleep 30
done
