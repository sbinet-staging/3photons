with Text_IO;

procedure Unit is
   use Text_IO;
   学生数 : Integer := 0;

begin
   Put_Line (Integer'Image (学生数));
   学生数 := 学生数 + 1;
   Put_Line  (Integer'Image (学生数));
end Unit;
