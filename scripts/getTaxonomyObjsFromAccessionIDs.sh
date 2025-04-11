# a file of accession ids, separated by new line characters
ACCESSION_ID_PATH=$1
# the output file
TAXONOMY_OBJS_PATH=$2
INTERMEDIATE_PATH="/tmp/tree.txt"
if [ -z "$ACCESSION_ID_PATH" ] || [ -z "$INTERMEDIATE_PATH" ]; then
  echo "Incorrect number of arguments. Need at least two paths:"
  echo -e "\tInput Path"
  echo -e "\tOutput Path"
  exit 1
fi
# clear the contents of the output file
truncate -s 0 $INTERMEDIATE_PATH
# append a newline character to the end of the file if it doesn't end with one
if [[ $(tail -c 1 $ACCESSION_ID_PATH) != "" ]]; then
  echo "" >> "$ACCESSION_ID_PATH"
fi
function processAccessionIDs() { echo $numIDs
  # remove the trailing comma and space
  ACCESSION_IDS=${ACCESSION_IDS%, }
  # get the taxonomy objects associated with each accession ID
  TAXONOMY_OBJS=$( ncbi-taxonomist map -edb assembly --accessions $ACCESSION_IDS | ncbi-taxonomist resolve -m )
  # print the taxonomyObjs to the output file
  echo -e "$TAXONOMY_OBJS" >> $INTERMEDIATE_PATH
}
# loop through each accession id and append to a variable
# process the IDs in large chunks
numIDs=0
ACCESSION_IDS=""
while read ACCESSION_ID; do
  ACCESSION_IDS+="$ACCESSION_ID, "
  numIDs=$((numIDs + 1))
  if [ $(( $numIDs % 300 )) -eq 0 ] ; then
    processAccessionIDs $numIDs $ACCESSION_IDS
    ACCESSION_IDS=""
  fi

#done < "$ACCESSION_ID_PATH"
processAccessionIDs $ACCESSION_IDS

mv "$TAXONOMY_OBJS_PATH" "$INTERMEDIATE_PATH"

echo "done."
