for ipy in  $@
do
    jupyter nbconvert --to python $ipy
done
