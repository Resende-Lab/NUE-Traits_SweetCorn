#Replace oldest hybrid parent with parent of best hybrid from YT6

Elite_private = Parents_private_update #c(Elite_private,Parents_private_update)
Elite_private = selectInd(Elite_private, nInd = nElite, selectTop = TRUE, use="gv",trait=selIndex, b=b)


#Update testers
Tester1_private = Elite_private[1:nTester1]
Tester2_private = Elite_private[1:nTester2]
Tester3_private = Elite_private[1:nTester3]
