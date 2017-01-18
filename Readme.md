# Read me file

## Smaller header   
### Even smaller header   

# Annotating the ecological genomics course 

Create a list:   
* first   
* second   
* third   

Syntax for creating a list   

```
* Item 1   
* Item 2   
* Item 3   
```

Example of embedding a URL:   
Syntax:   
```
[hyperlinked words](http://www.pinkbike.com)
```
Implementation:

This is a cool [website](http://www.pinkbike.com).

Example of embedding an image

Syntax:

```
![](URL)
```

This is a picture.   
![](https://cloud.githubusercontent.com/assets/21958390/22071930/b9347194-dd6e-11e6-987a-0735adb739c4.jpeg)

#packages for reading in data 
library(data.table) 

dat<-fread("csvfile")
#making a table
knit::kable(dat)
