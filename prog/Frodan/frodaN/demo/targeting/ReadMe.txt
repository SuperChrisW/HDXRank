11. Set environment variable to point to the directory where frodaN package is installed. For example: 
export FRODANHOME=/home/user/frodaN/

2. Run prepossesing script:
python $FRODANHOME/preprocess.py -i 1CFC.pdb -t 1CFD.pdb -o options.xml

3. Run frodaN: 
$FRODANHOME/bin/main options.xml

