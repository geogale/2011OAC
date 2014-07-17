#####    #####     #####    #####     #####    #####     #####    #####     #####    #####     #####
# Code to Create the 2011 Area Classification for Output Areas (2011 OAC)
# Feel free to share and reuse with attribution
# Ancillary 5 - Combined CSV files in Microsoft Excel
#####    #####     #####    #####     #####    #####     #####    #####     #####    #####     #####

####################################################################################################
# Import CSV #######################################################################################
####################################################################################################

Option Explicit
Sub ImportCSVs()
'Summary:   Import all CSV files from a folder into separate sheets
'           named for the CSV filenames
Dim fPath   As String
Dim fCSV    As String
Dim wbCSV   As Workbook
Dim wbMST   As Workbook

Set wbMST = ThisWorkbook
fPath = "C:\Users\Chris\Documents\2011 OAC\FOLDERNAME\"                  'path to CSV files, include the final \
Application.ScreenUpdating = False  'speed up macro
Application.DisplayAlerts = False   'no error messages, take default answers
fCSV = Dir(fPath & "*.csv")         'start the CSV file listing

    On Error Resume Next
    Do While Len(fCSV) > 0
        Set wbCSV = Workbooks.Open(fPath & fCSV)                    'open a CSV file
        wbMST.Sheets(ActiveSheet.Name).Delete                       'delete sheet if it exists
        ActiveSheet.Move After:=wbMST.Sheets(wbMST.Sheets.Count)    'move new sheet into Mstr
        Columns.Autofit             'clean up display 
        fCSV = Dir                  'ready next CSV
    Loop
  
Application.ScreenUpdating = True
Set wbCSV = Nothing
End Sub

####################################################################################################
# Freeze Top Row ###################################################################################
####################################################################################################

Sub FreezeTopRow()
    Application.ScreenUpdating = False
    For Each ws In ThisWorkbook.Worksheets
        ws.Activate
        With ActiveWindow
            .SplitColumn = 0: .SplitRow = 1
            .FreezePanes = True
        End With
    Next
End Sub

####################################################################################################
# Justification ####################################################################################
####################################################################################################

Sub Justification()
    Application.ScreenUpdating = False
    For Each ws In ThisWorkbook.Worksheets
        ws.Activate
        With ActiveWindow
    Cells.Select
    With Selection.Font
        .Name = "Cambria"
        .Size = 11
        .Strikethrough = False
        .Superscript = False
        .Subscript = False
        .OutlineFont = False
        .Shadow = False
        .Underline = xlUnderlineStyleNone
        .ThemeColor = xlThemeColorLight1
        .TintAndShade = 0
        .ThemeFont = xlThemeFontNone
    End With
    Rows("1:1").Select
    Selection.Font.Bold = True
    Columns("A:P").Select
    Columns("A:P").EntireColumn.AutoFit
    ActiveWindow.ScrollColumn = 11
    ActiveWindow.ScrollColumn = 10
    ActiveWindow.ScrollColumn = 9
    ActiveWindow.ScrollColumn = 8
    ActiveWindow.ScrollColumn = 7
    ActiveWindow.ScrollColumn = 6
    ActiveWindow.ScrollColumn = 5
    ActiveWindow.ScrollColumn = 4
    ActiveWindow.ScrollColumn = 1
    Columns("A:A").Select
    With Selection
        .HorizontalAlignment = xlRight
        .VerticalAlignment = xlBottom
        .WrapText = False
        .Orientation = 0
        .AddIndent = False
        .IndentLevel = 0
        .ShrinkToFit = False
        .ReadingOrder = xlContext
        .MergeCells = False
    End With
    Columns("B:B").Select
    With Selection
        .HorizontalAlignment = xlLeft
        .VerticalAlignment = xlBottom
        .WrapText = False
        .Orientation = 0
        .AddIndent = False
        .IndentLevel = 0
        .ShrinkToFit = False
        .ReadingOrder = xlContext
        .MergeCells = False
    End With
    Columns("C:C").Select
    With Selection
        .HorizontalAlignment = xlRight
        .VerticalAlignment = xlBottom
        .WrapText = False
        .Orientation = 0
        .AddIndent = False
        .IndentLevel = 0
        .ShrinkToFit = False
        .ReadingOrder = xlContext
        .MergeCells = False
    End With
    Columns("D:D").Select
    With Selection
        .HorizontalAlignment = xlLeft
        .VerticalAlignment = xlBottom
        .WrapText = False
        .Orientation = 0
        .AddIndent = False
        .IndentLevel = 0
        .ShrinkToFit = False
        .ReadingOrder = xlContext
        .MergeCells = False
    End With
    Columns("E:E").Select
    With Selection
        .HorizontalAlignment = xlRight
        .VerticalAlignment = xlBottom
        .WrapText = False
        .Orientation = 0
        .AddIndent = False
        .IndentLevel = 0
        .ShrinkToFit = False
        .ReadingOrder = xlContext
        .MergeCells = False
    End With
    ActiveWindow.SmallScroll ToRight:=4
    Columns("F:F").Select
    With Selection
        .HorizontalAlignment = xlLeft
        .VerticalAlignment = xlBottom
        .WrapText = False
        .Orientation = 0
        .AddIndent = False
        .IndentLevel = 0
        .ShrinkToFit = False
        .ReadingOrder = xlContext
        .MergeCells = False
    End With
    Columns("G:G").Select
    With Selection
        .HorizontalAlignment = xlRight
        .VerticalAlignment = xlBottom
        .WrapText = False
        .Orientation = 0
        .AddIndent = False
        .IndentLevel = 0
        .ShrinkToFit = False
        .ReadingOrder = xlContext
        .MergeCells = False
    End With
    Columns("H:H").Select
    With Selection
        .HorizontalAlignment = xlLeft
        .VerticalAlignment = xlBottom
        .WrapText = False
        .Orientation = 0
        .AddIndent = False
        .IndentLevel = 0
        .ShrinkToFit = False
        .ReadingOrder = xlContext
        .MergeCells = False
    End With
    ActiveWindow.SmallScroll ToRight:=2
    Columns("I:I").Select
    With Selection
        .HorizontalAlignment = xlRight
        .VerticalAlignment = xlBottom
        .WrapText = False
        .Orientation = 0
        .AddIndent = False
        .IndentLevel = 0
        .ShrinkToFit = False
        .ReadingOrder = xlContext
        .MergeCells = False
    End With
    Columns("J:J").Select
    With Selection
        .HorizontalAlignment = xlLeft
        .VerticalAlignment = xlBottom
        .WrapText = False
        .Orientation = 0
        .AddIndent = False
        .IndentLevel = 0
        .ShrinkToFit = False
        .ReadingOrder = xlContext
        .MergeCells = False
    End With
    ActiveWindow.SmallScroll ToRight:=3
    Columns("K:K").Select
    With Selection
        .HorizontalAlignment = xlRight
        .VerticalAlignment = xlBottom
        .WrapText = False
        .Orientation = 0
        .AddIndent = False
        .IndentLevel = 0
        .ShrinkToFit = False
        .ReadingOrder = xlContext
        .MergeCells = False
    End With
    Columns("L:L").Select
    With Selection
        .HorizontalAlignment = xlLeft
        .VerticalAlignment = xlBottom
        .WrapText = False
        .Orientation = 0
        .AddIndent = False
        .IndentLevel = 0
        .ShrinkToFit = False
        .ReadingOrder = xlContext
        .MergeCells = False
    End With
    Columns("M:M").Select
    With Selection
        .HorizontalAlignment = xlRight
        .VerticalAlignment = xlBottom
        .WrapText = False
        .Orientation = 0
        .AddIndent = False
        .IndentLevel = 0
        .ShrinkToFit = False
        .ReadingOrder = xlContext
        .MergeCells = False
    End With
    Columns("N:N").Select
    With Selection
        .HorizontalAlignment = xlLeft
        .VerticalAlignment = xlBottom
        .WrapText = False
        .Orientation = 0
        .AddIndent = False
        .IndentLevel = 0
        .ShrinkToFit = False
        .ReadingOrder = xlContext
        .MergeCells = False
    End With
    ActiveWindow.SmallScroll ToRight:=4
    Columns("O:O").Select
    With Selection
        .HorizontalAlignment = xlRight
        .VerticalAlignment = xlBottom
        .WrapText = False
        .Orientation = 0
        .AddIndent = False
        .IndentLevel = 0
        .ShrinkToFit = False
        .ReadingOrder = xlContext
        .MergeCells = False
    End With
    Columns("P:P").Select
    With Selection
        .HorizontalAlignment = xlLeft
        .VerticalAlignment = xlBottom
        .WrapText = False
        .Orientation = 0
        .AddIndent = False
        .IndentLevel = 0
        .ShrinkToFit = False
        .ReadingOrder = xlContext
        .MergeCells = False
    End With
    ActiveWindow.ScrollColumn = 13
    ActiveWindow.ScrollColumn = 12
    ActiveWindow.ScrollColumn = 9
    ActiveWindow.ScrollColumn = 8
    ActiveWindow.ScrollColumn = 7
    ActiveWindow.ScrollColumn = 6
    ActiveWindow.ScrollColumn = 5
    ActiveWindow.ScrollColumn = 4
    ActiveWindow.ScrollColumn = 3
    ActiveWindow.ScrollColumn = 2
    ActiveWindow.ScrollColumn = 1
    Range("A1").Select
        End With
    Next
End Sub

####################################################################################################
# Protect All Sheets ###############################################################################
####################################################################################################

Sub ProtectAll()
	Dim wSheet As Worksheet
	For Each wSheet In Worksheets
		wSheet.Protect DrawingObjects:=True, Contents:=True, Scenarios:=True _
    	    , AllowFormattingCells:=True, AllowFormattingColumns:=True, _
     	   AllowFormattingRows:=True, AllowSorting:=True, AllowFiltering:=True
	Next wSheet
End Sub
