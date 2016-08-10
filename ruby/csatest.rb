require 'csa'

csa = CSA.new(ARGV[0],ARGV[1])
for i in 0..9
  j = csa.lookup(i)
  print("SA[", i, "] = ", j, "\n")
end

while 1 do
  print("input key ")
  key = STDIN.gets().chomp
  puts "key => #{key.inspect}"
  break if key.empty?

  range = csa.search(key)
  if range != nil then
    puts "[#{range[0]}, #{range[1]}]"
    for i in range[0]..range[1]
      j = csa.lookup(i)
      puts "SA[#{i}]=#{j}: #{csa.text([j-20,j+20])}"
#      puts "#{j}: #{csa.substring(i, 20)}"
    end
  else
    puts "not found"
  end
end
