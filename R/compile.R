## LEGO
system("gcc -o LEGO_pre LEGO_pre.c -lm")
system("gcc -oLEGO LEGO.c -lm")
system("gcc -o LEGO_mul LEGO_mul.c -lm")
system("gcc -o LEGO_bg LEGO_bg.c -lm")
system("gcc -o LEGO_mul_bg LEGO_mul_bg.c -lm")
system("g++ -o extract_info extract_info.cpp -lm")
system("g++ -o extract_info_gs extract_info_gs.cpp -lm")
##
inp = "iNP.c"
new_code = paste(inp, ".new", sep="")
new_code = inp

system("gcc -o iNP iNP.c -lm")

code1 = "MSG.cpp"
code2 = "VM.cpp"
new_code = paste(code1, ".new", sep="")
new_code <- code1
system("g++ code1 -o MSG.out") 
system("g++ code2 -o VM.out")
system("mv MSG.out") 
system("mv VM.out")

