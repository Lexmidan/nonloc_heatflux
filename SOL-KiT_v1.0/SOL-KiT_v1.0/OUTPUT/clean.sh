#Cleans all .txt files from output folders
echo 'Cleaning output...' 
for dir in */
do
    rm -f ./${dir}*.txt
done
echo 'Removing Status file if present'
rm -f ./STATUS.txt
echo 'Done'

