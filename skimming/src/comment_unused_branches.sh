#! /bin/zsh

SELECTOR_H="GeneralSkimmer.h"
SELECTOR_C="GeneralSkimmer.cpp"

# get list of variable names
TREE_VARS=$(grep "TBranch        \*" $SELECTOR_H | cut -c22- | cut -d';' -f1);
# convert to an array
saveIFS=$IFS; IFS=$'\n'; array=($(echo $TREE_VARS)); IFS=$saveIFS;

# check if variable is in selector implementation and comment if not
for var in $array; 
do
  if $(grep -q $var $SELECTOR_C)
  then
    echo $var" is used"
  else
    echo $var" not found"
    pattern="/"$var"/ s?^?//?"
    echo $pattern
    # comment if var in line
    sed -e $pattern -i $SELECTOR_H
    # fix multiple comment (names contained in others)
    sed -e "s?^////?//?" -i $SELECTOR_H
  fi
done
