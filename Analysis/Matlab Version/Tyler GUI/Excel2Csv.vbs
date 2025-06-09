if WScript.Arguments.Count < 5 Then
    WScript.Echo "Please specify the the sourcefile, sheetname1, destinationfile1, sheetname2, destinationfile2"
    Wscript.Quit
End If

csv_format = 6

Set objFSO = CreateObject("Scripting.FileSystemObject")

src_file = objFSO.GetAbsolutePathName(Wscript.Arguments.Item(0))
dest_file1 = objFSO.GetAbsolutePathName(WScript.Arguments.Item(2))
dest_file2 = objFSO.GetAbsolutePathName(WScript.Arguments.Item(4))

Dim oExcel
Set oExcel = CreateObject("Excel.Application")

Dim oBook
Set oBook = oExcel.Workbooks.Open(src_file)

oBook.Sheets(WScript.Arguments.Item(1)).Select
oBook.SaveAs dest_file1, csv_format

oBook.Sheets(WScript.Arguments.Item(3)).Select
oBook.SaveAs dest_file2, csv_format

oBook.Close False
oExcel.Quit