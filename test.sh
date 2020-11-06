while IFS= read -r line; do

    echo "$line"
    cp submission/"$line/"*.cpp src/
    # cp submission/"$line/"*.h include/

    mkdir build"_$line"
    cd build"_$line"
    cmake .. -DCMAKE_BUILD_TYPE=Release
    make -j8

    cd ../

    # exit 0

done < hw2-namelist.txt


# for d in * ; do

#     # echo "$d"
#     # # # # cp base_include/* include/
#     # cp $d/*.cpp src/
#     # # # cp submission/"$d/"*.h include/
#     # # mkdir ../build"_$d"
#     # # cd ../build"_$d"
#     # # cmake .. -DCMAKE_BUILD_TYPE=Release
#     # # make -j8
#     # # cd ../submission

#     # echo $d/

#     # # cd $d/

#     # name = basename $d
#     echo $d >> ../hw2-namelist.txt

#     # exit 0

# done