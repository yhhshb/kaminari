function debian_realpath() {
    f=$@;
    if [ -d "$f" ]; then
        base="";
        dir="$f";
    else
        base="/$(basename "$f")";
        dir=$(dirname "$f");
    fi
    dir=$(cd "$dir" && /bin/pwd);
    echo "$dir$base";
}

script_name=$0
script_full_path=$(dirname $(debian_realpath "$0"))
lphash_full_path=$(dirname "$script_full_path")

salmonella10_folder="$lphash_full_path/data/salmonella10"

mkdir -p $salmonella10_folder

SAL_AA7743AA=SAL_AA7743AA.fasta.gz
SAL_BA0010AA=SAL_BA0010AA.fasta.gz
SAL_CA3280AA=SAL_CA3280AA.fasta.gz
SAL_FA0063AA=SAL_FA0063AA.fasta.gz
SAL_FA6579AA=SAL_FA6579AA.fasta.gz
SAL_GA5038AA=SAL_GA5038AA.fasta.gz
SAL_HA1487AA=SAL_HA1487AA.fasta.gz
SAL_HA3099AA=SAL_HA3099AA.fasta.gz
SAL_HA8439AA=SAL_HA8439AA.fasta.gz
SAL_HA8462AA=SAL_HA8462AA.fasta.gz

wget "https://github.com/jermp/fulgor/raw/main/test_data/salmonella_10/$SAL_AA7743AA"
wget "https://github.com/jermp/fulgor/raw/main/test_data/salmonella_10/$SAL_BA0010AA"
wget "https://github.com/jermp/fulgor/raw/main/test_data/salmonella_10/$SAL_CA3280AA"
wget "https://github.com/jermp/fulgor/raw/main/test_data/salmonella_10/$SAL_FA0063AA"
wget "https://github.com/jermp/fulgor/raw/main/test_data/salmonella_10/$SAL_FA6579AA"
wget "https://github.com/jermp/fulgor/raw/main/test_data/salmonella_10/$SAL_GA5038AA"
wget "https://github.com/jermp/fulgor/raw/main/test_data/salmonella_10/$SAL_HA1487AA"
wget "https://github.com/jermp/fulgor/raw/main/test_data/salmonella_10/$SAL_HA3099AA"
wget "https://github.com/jermp/fulgor/raw/main/test_data/salmonella_10/$SAL_HA8439AA"
wget "https://github.com/jermp/fulgor/raw/main/test_data/salmonella_10/$SAL_HA8462AA"

mv $SAL_AA7743AA $salmonella10_folder/$SAL_AA7743AA
mv $SAL_BA0010AA $salmonella10_folder/$SAL_BA0010AA
mv $SAL_CA3280AA $salmonella10_folder/$SAL_CA3280AA
mv $SAL_FA0063AA $salmonella10_folder/$SAL_FA0063AA
mv $SAL_FA6579AA $salmonella10_folder/$SAL_FA6579AA
mv $SAL_GA5038AA $salmonella10_folder/$SAL_GA5038AA
mv $SAL_HA1487AA $salmonella10_folder/$SAL_HA1487AA
mv $SAL_HA3099AA $salmonella10_folder/$SAL_HA3099AA
mv $SAL_HA8439AA $salmonella10_folder/$SAL_HA8439AA
mv $SAL_HA8462AA $salmonella10_folder/$SAL_HA8462AA
