
for i in `seq 0 214`;do
  if [ "$i" -lt "109" ]; then
    sed -e "s/nframes/1001/g" -e "s/num/$i/g" unwraptemp.in > unwrap.in
  else
    sed -e "s/nframes/2001/g" -e "s/num/$i/g" unwraptemp.in > unwrap.in
  fi
  ../unwrap unwrap.in
done

