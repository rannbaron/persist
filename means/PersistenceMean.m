function [ mean ] = PersistenceMean( cell_diagrams, epsilon, zigzag )

    mean = kPersistenceMean( cell_diagrams, epsilon, zigzag, 1);

end