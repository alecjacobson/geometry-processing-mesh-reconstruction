if [ $# -le 0 ]; then
  echo "Invalid argument, expect: ./grade.sh [id]"
  exit 1
fi

cd build"_$1"
./mesh-reconstruction
./mesh-reconstruction ../data/elephant.pwn

# # cp base_include/* include/
# cp submission/"$1/"*.cpp src/
# cp submission/"$1/"*.h include/
# mkdir build"_$1"
# cd build"_$1"
# cmake .. -DCMAKE_BUILD_TYPE=Release
# make -j8


# ./a1-mass-spring-1d
# ./a1-mass-spring-1d rk
# ./a1-mass-spring-1d be
# ./a1-mass-spring-1d se
# cp ./a1-mass-spring-1d ./a1-mass-spring-1d"_$1"
# cd ..
