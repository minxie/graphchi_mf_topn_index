export GRAPHCHI_ROOT=../
for m in '50'
do
	for l in  '1e-3'
  	do
		../bin/myapps/sgd --D=$m --training=../../data/ml-10M100K/ratings.dat_mm --result=../../result/gc_movielens_res \
		  --sgd_gamma=1e-4 --sgd_lambda=1e-4 --max_iter=100 --n_top=50 --quiet=1 --tol=1e-3
    	done
done
