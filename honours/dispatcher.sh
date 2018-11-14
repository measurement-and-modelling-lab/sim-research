sed 1d "conditions.csv" | while read line
do
arguments=${line//\,/ }          ## Replace commas with spaces
./a.out $arguments >> output.csv ## Feed arguments to C++ and redirect output to file
done
