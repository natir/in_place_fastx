var searchIndex = JSON.parse('{\
"in_place_fastx":{"doc":"","t":[17,0,0,4,13,13,13,13,13,13,6,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,12,12,12,3,3,0,11,11,11,11,12,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,0,12,12,12,11,11,11,11,11,11,3,3,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,0,0,8,11,11,10,11,8,11,11],"n":["DEFAULT_BLOCKSIZE","error","fastq","Error","MapFile","MetaDataFile","NoNewLineInBlock","NotAFastqFile","OpenFile","PartialRecord","Result","borrow","borrow_mut","deref","deref_mut","drop","fmt","fmt","from","init","into","source","to_string","try_from","try_into","type_id","source","source","source","Block","Record","block","borrow","borrow","borrow_mut","borrow_mut","comment","data","deref","deref","deref_mut","deref_mut","drop","drop","from","from","init","init","into","into","is_empty","len","new","parser","plus","quality","sequence","try_from","try_from","try_into","try_into","type_id","type_id","Producer","Reader","borrow","borrow","borrow_mut","borrow_mut","deref","deref","deref_mut","deref_mut","drop","drop","from","from","init","init","into","into","into_iter","new","new","next","next_block","next_record","try_from","try_from","try_into","try_into","type_id","type_id","with_blocksize","sequential","shared_state","Sequential","block","parse","record","with_blocksize","SharedState","parse","with_blocksize"],"q":["in_place_fastx","","","in_place_fastx::error","","","","","","","","","","","","","","","","","","","","","","","in_place_fastx::error::Error","","","in_place_fastx::fastq","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","in_place_fastx::fastq::block","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","in_place_fastx::fastq::parser","","in_place_fastx::fastq::parser::sequential","","","","","in_place_fastx::fastq::parser::shared_state","",""],"d":["","","","","","","","","","","","","","","","","","","","","","","","","","","","","","Block reperesent a section of file memory mapped in file","Record store a fastq record in public field","Struct that extract part of file (called block) and read …","","","","","","Acces to data owned by block","","","","","","","","","","","","","Return true if the block is empty","Get length of block","Create a new Block","","","","","","","","","","","Struct that produce a Block of file, this block contains …","Struct that read Block and produce Record","","","","","","","","","","","","","","","","","","Build a Block producer, default Block size is 2^16 bytes.","Create a new Block reader from Block get in parameter","","Get the next Block, all Block contains almost one record","Produce Record until block is empty","","","","","","","Build a Block producer, with a specific Block size …","","Struct that extract part of file (called block), each …","Trait allow sequential parsing of fastq","Method call to parse a block","Parse file indicate by path with default blocksize […","Method call to parse a record","Parse file indicate by path with selected blocksize","Trait allow parallel parsing of block.","Parse file indicate by path with default blocksize […","Parse file indicate by path with selected blocksize"],"i":[0,0,0,0,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,3,4,0,0,0,5,6,5,6,5,6,5,6,5,6,5,6,5,6,5,6,5,6,6,6,6,0,5,5,5,5,6,5,6,5,6,0,0,7,8,7,8,7,8,7,8,7,8,7,8,7,8,7,8,7,7,8,7,7,8,7,8,7,8,7,8,7,0,0,0,9,9,9,9,0,10,10],"f":[null,null,null,null,null,null,null,null,null,null,null,[[]],[[]],[[["usize",15]]],[[["usize",15]]],[[["usize",15]]],[[["formatter",3]],["result",6]],[[["formatter",3]],["result",6]],[[]],[[],["usize",15]],[[]],[[],[["error",8],["option",4]]],[[],["string",3]],[[],["result",4]],[[],["result",4]],[[],["typeid",3]],null,null,null,null,null,null,[[]],[[]],[[]],[[]],null,[[]],[[["usize",15]]],[[["usize",15]]],[[["usize",15]]],[[["usize",15]]],[[["usize",15]]],[[["usize",15]]],[[]],[[]],[[],["usize",15]],[[],["usize",15]],[[]],[[]],[[],["bool",15]],[[],["usize",15]],[[["usize",15],["mmap",3]]],null,null,null,null,[[],["result",4]],[[],["result",4]],[[],["result",4]],[[],["result",4]],[[],["typeid",3]],[[],["typeid",3]],null,null,[[]],[[]],[[]],[[]],[[["usize",15]]],[[["usize",15]]],[[["usize",15]]],[[["usize",15]]],[[["usize",15]]],[[["usize",15]]],[[]],[[]],[[],["usize",15]],[[],["usize",15]],[[]],[[]],[[]],[[],["result",6]],[[["block",3]]],[[],["option",4]],[[],[["result",6],["option",4]]],[[],[["option",4],["result",6]]],[[],["result",4]],[[],["result",4]],[[],["result",4]],[[],["result",4]],[[],["typeid",3]],[[],["typeid",3]],[[["u64",15]],[["error",4],["result",4]]],null,null,null,[[["block",3]],["result",6]],[[],["result",6]],[[["record",3]]],[[["u64",15]],["result",6]],null,[[],["result",6]],[[["u64",15]],["result",6]]],"p":[[4,"Error"],[13,"MetaDataFile"],[13,"OpenFile"],[13,"MapFile"],[3,"Record"],[3,"Block"],[3,"Producer"],[3,"Reader"],[8,"Sequential"],[8,"SharedState"]]}\
}');
if (window.initSearch) {window.initSearch(searchIndex)};