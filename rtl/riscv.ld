/* RISC-V 链接脚本 - 用于 AV2 SoC */

MEMORY
{
    ROM (rx)  : ORIGIN = 0x00000000, LENGTH = 64K
    RAM (rwx) : ORIGIN = 0x20000000, LENGTH = 64K
}

SECTIONS
{
    /* 程序代码放在 ROM */
    .text :
    {
        *(.text)
        *(.text.*)
        *(.rodata)
        *(.rodata.*)
    } > ROM

    /* 变量和堆栈放在 RAM */
    .data :
    {
        *(.data)
        *(.data.*)
        *(.sdata)
        *(.sdata.*)
    } > RAM AT > ROM

    .bss :
    {
        *(.bss)
        *(.bss.*)
        *(.sbss)
        *(.sbss.*)
        *(COMMON)
    } > RAM

    /* 定义堆栈顶 (在 RAM 末尾) */
    _stack_top = ORIGIN(RAM) + LENGTH(RAM);
}
