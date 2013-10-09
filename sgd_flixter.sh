export GRAPHCHI_ROOT=../
for m in '50'
do
	for l in  '1e-3'
  	do
		../bin/myapps/sgd --D=$m --training=../../data/flixter/ratings_flixter_mm \
		  --sgd_gamma=1e-4 --sgd_lambda=1e-4 --max_iter=6 --n_top=50 --quiet=1
    	done
done
