
/* WARNING: Globals starting with '_' overlap smaller symbols at the same address */

undefined8 main(int param_1,undefined8 *param_2)

{
  char *pcVar1;
  ulong uVar2;
  undefined8 uVar3;
  FILE *__stream;
  char *__ptr;
  long lVar4;
  long lVar5;
  size_t __n;
  long in_FS_OFFSET;
  ulong local_70;
  undefined1 local_68 [16];
  undefined1 local_58 [16];
  char *local_48;
  char *local_40;
  char *local_38;
  long local_30;
  
  local_30 = *(long *)(in_FS_OFFSET + 0x28);
  if (param_1 == 2) {
    __stream = fopen((char *)param_2[1],"rb");
    if (__stream == (FILE *)0x0) {
      perror("Failed to open file");
      goto LAB_0010015e;
    }
    fread(&local_70,8,1,__stream);
    uVar2 = local_70;
    if ((long)local_70 < 0) {
      if (local_30 == *(long *)(in_FS_OFFSET + 0x28)) {
        std::__throw_length_error("cannot create std::vector larger than max_size()");
      }
      goto LAB_00100308;
    }
    if (local_70 == 0) {
      __ptr = (char *)0x0;
      local_48 = (char *)0x0;
      local_38 = (char *)0x0;
      local_40 = (char *)0x0;
    }
    else {
      __ptr = (char *)operator_new(local_70);
      pcVar1 = __ptr + uVar2;
      *__ptr = '\0';
      __n = uVar2 - 1;
      local_40 = __ptr + 1;
      local_48 = __ptr;
      local_38 = pcVar1;
      if (__n != 0) {
        memset(__ptr + 1,0,__n);
        local_40 = pcVar1;
      }
    }
                    /* try { // try from 00100210 to 001002b2 has its CatchHandler @ 0010030d */
    fread(__ptr,1,local_70,__stream);
    fclose(__stream);
    local_68 = (undefined1  [16])0x0;
    local_58 = (undefined1  [16])0x0;
    lVar4 = std::chrono::_V2::system_clock::now();
    simulate(local_70,local_48,(complex *)local_68,(complex *)local_58);
    lVar5 = std::chrono::_V2::system_clock::now();
    __printf_chk(local_68._0_8_,local_68._8_8_,local_58._0_8_,local_58._8_8_,2,
                 "Final state: alpha = %.12f + %.12fi, beta = %.12f + %.12fi\n");
    __printf_chk((double)(lVar5 - lVar4) / 1000000.0,2,"Time taken: %.2f ms\n");
    std::vector<>::~vector((vector<> *)&local_48);
    uVar3 = 0;
  }
  else {
    __fprintf_chk(_stderr,2,"Usage: %s <input_file>\n",*param_2);
LAB_0010015e:
    uVar3 = 1;
  }
  if (local_30 == *(long *)(in_FS_OFFSET + 0x28)) {
    return uVar3;
  }
LAB_00100308:
                    /* WARNING: Subroutine does not return */
  __stack_chk_fail();
}

