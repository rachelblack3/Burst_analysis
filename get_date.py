def get_date_string(single_date):

    date_string= str(single_date.strftime("%Y%m%d"))

    if (single_date.day <10):
        day = "0"+str(single_date.day)
    else:
        day = str(single_date.day)


    if (single_date.month<10):
        month = "0"+str(single_date.month)
    else:
        month = str(single_date.month)
    

    year = str(single_date.year)
    
    return date_string,year,month,day

   
