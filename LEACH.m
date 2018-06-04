clear;

%Alan boyutlar� x ve y metre olarak
xm = 100;
ym = 100;

%Baz istasyonun x ve y kordinatlar�
sink.x = 1.5 * xm;
sink.y = 0.5 * ym;

%Alandaki node say�s�
n = 200;

%Bir nodun clusterhead se�ilmesi i�in gereken olas�l�k
p = 0.2;

%Ba�lang�� enerjileri Joule cinsinden
Eo = 0.1; %Her bir nodun ba�lang�� enerjisi
%Eelec=Etx=Erx Veri iletmek veya almak i�in elektronik cihaz taraf�ndan harcanan enerji
ETX = 50 * 0.000000001; %Veri iletim enerjisi
ERX = 50 * 0.000000001; %Veri alma enerjisi

Efs = 10 * 0.000000000001; %Da��t�m s�ras�nda kaybolan enerji enerjisi
Emp = 0.0013 * 0.000000000001; %Enerji transferinin amplifikasyon katsay�s�
EDA = 5 * 0.000000001; %CH taraf�ndan veri toplama s�ras�nda harcad�klar� enerjisi

k = 4000; %Veri paket boyutu

%Advanced node y�zdesi(y�ksek enerjili node y�zdesi)
m = 0.0;
%\alpha - advanced nodelara enerji atarken kullan�lan kat say�
a = 1;

%Maksimum iterasyon say�s�
rmax = 100;

%BS'den al�c� d���me olan uzakl�k, Emp veya Efs'e ba�l�d�r.
do = sqrt(Efs / Emp);

%Sensor network�n olu�turulmas�
figure(1);
hold off;
for i = 1:1:n
    S(i).xd = rand(1, 1) * xm; %i ninci nodun x kordinat�n�  xm degeri aral�g�nda �retip (S node dizisi)
    XR(i) = S(i).xd; %XR dizisine x kordinat�n� atad�k
    S(i).yd = rand(1, 1) * ym; %i ninci nodun y kordinat�n�  yd degeri aral�g�nda �retip (S node dizisi)
    YR(i) = S(i).yd; %YR dizisine y kordinat�n� atad�k
    S(i).G = 0; %Clusterhead (CH) olabilecek advencad(enerjisi en y�ksek olanlar) node lar�n yada normal nodelar�n k�mesi
    S(i).type = 'N'; %Ba�lang��ta hi� cluster head olmad� i�in nodelar�n tipini n olarak i�aretledik yani normal node
 
    %Belirtilen advance node oran�na g�re, Normal node lar�n se�ilmesi ve enerjilerinin atanmas�
    if (i >= m * n + 1)
        S(i).E = Eo;
        S(i).IS_ADVANCED = 0; %Advanced olmayacak
        plot(S(i).xd, S(i).yd, 'o'); %Normal node u ekrana �izdirdik
        hold on;
    end
    %Belirtilen advance node oran�na g�re ilk x tanesini advence node olarak se� ve enerjilerinin atanmas�
    if (i < m * n + 1)
        S(i).E = Eo * (1 + a)
        S(i).IS_ADVANCED = 1; %Advanced olacak
        plot(S(i).xd, S(i).yd, '+'); %Y�ksek enerjili(advanced) nodu ekrana �izdirdik
        hold on;
    end
end

%S dizinde node listemiz var node dizisinin en sonuna  baz istesyonunun kordinatlar�n� ekledik
S(n + 1).xd = sink.x;
S(n + 1).yd = sink.y;
plot(S(n + 1).xd, S(n + 1).yd, 'x'); %Baz istasyonunu ekrana x olarak �izdirdik

%�lk iterasyon i�in tan�mlar

%Cluster head sayac�
countCHs = 0;
%Round daki cluster head say�s�n� tutan degi�ken
rcountCHs = 0;

%�lk �l�m�n ger�ekle�ip ger�ekle�medi�inden haberdar oldugumuz degi�ken
flag_first_dead = 0;

for r = 0:1:rmax
    %r
 
    %Operation for epoch
    if (mod(r, round(1 / p)) == 0)
        for i = 1:1:n
            S(i).G = 0;
            S(i).cl = 0;
        end
    end
 
    hold off;
 
    %Her round ba�lad���nda toplam �l� node say�s� tutuyor, her round ba�lad���nda s�f�rlan�yor
    total_dead = 0;
    %Her round ba�lad���nda Advanced �l� node say�s�, her round ba�lad���nda s�f�rlan�yor
    dead_a = 0;
    %Her round ba�lad���nda Normal �l� node say�s�, her round ba�lad���nda s�f�rlan�yor
    dead_n = 0;
 
    %BS ve CH ye g�nderilen paket sayac�
    packets_TO_BS = 0;
    packets_TO_CH = 0;
    %Roundlarda CH ve BS a g�nderilen paket sayac� dizi
    PACKETS_TO_CH(r + 1) = 0;
    PACKETS_TO_BS(r + 1) = 0;
 
    figure(1);
 
    for i = 1:1:n
        %�l� node var m� kontrol et ve say
        if (S(i).E <= 0)
            plot(S(i).xd, S(i).yd, 'red .'); %En son canl� kalanlar grafiginde �l� nodelar� k�rm�z� nokta olarak g�ster
            total_dead = total_dead + 1;
            if (S(i).IS_ADVANCED == 1)
                dead_a = dead_a + 1;
            end
            if (S(i).IS_ADVANCED == 0)
                dead_n = dead_n + 1;
            end
            hold on;
        end
     
        figure(2);
        %En son canl� kalan node lar
        if S(i).E > 0
            S(i).type = 'N';
            if (S(i).IS_ADVANCED == 0)
                plot(S(i).xd, S(i).yd, 'o');
            end
            if (S(i).IS_ADVANCED == 1)
                plot(S(i).xd, S(i).yd, '+');
            end
            hold on;
        end
    end
    plot(S(n + 1).xd, S(n + 1).yd, 'x'); %En son canl� kalanlar grafiginde baz istasyonunuda g�ster
 
    %STATISTICS dizisinin i�erisine toplam �l� advanced �l� normal �l� say�s�n� her bir round i�in atad�k
    STATISTICS(r + 1).DEAD = total_dead;
    STATISTICS(r + 1).DEAD_ADVANCED = dead_a;
    STATISTICS(r + 1).DEAD_NORMAL = dead_n;
 
    %DEAD dizisinin i�erisine toplam �l� advanced �l� normal �l� say�s�n� her bir round i�in atad�k
    DEAD(r + 1) = total_dead;
    DEAD_N(r + 1) = dead_n;
    DEAD_A(r + 1) = dead_a;
 
    %�lk nodun �ld��� roundu tespit ettik
    if (total_dead == 1)
        if (flag_first_dead == 0)
            first_dead = r;
            flag_first_dead = 1;
        end
    end
    %Cluster head sayac�, her round ba�lad�g�nda s�f�rlan�r
    countCHs = 0;
    cluster = 1; %Cluster say�s�n� 1 olarak at�yo ilk ba�ta
    for i = 1:1:n
        if (S(i).E > 0) %Elimizdeki node �l� degilse
            temp_rand = rand;
            if ((S(i).G) <= 0)
             
                %Clusterhead se�im
                if (temp_rand <= (p / (1 - p * mod(r, round(1 / p))))) %Clusterheadin se�im form�l�
                    countCHs = countCHs + 1; %Cluster say�s�n� 1 art�r
                    packets_TO_BS = packets_TO_BS + 1; %Baz istasyonuna g�nderilen paket say�s�n� 1 art�r
                    PACKETS_TO_BS(r + 1) = packets_TO_BS; %O round i�erisinde BS ye g�nderilen paket say�s�n� g�ncelliyor
                 
                    S(i).type = 'C'; %Se�ilen nodun tipini cluster olarak ayarla
                    S(i).G = round(1 / p) - 1; %Cluster head olarak se�ilen nodeun G de�erini ata
                    C(cluster).xd = S(i).xd; %Cluster�n x kordinat�n� C dizisine ata
                    C(cluster).yd = S(i).yd; %Cluster�n y kordinat�n� C dizisine ata
                    plot(S(i).xd, S(i).yd, 'k*'); %En son yani figure 10 a clusterheadlari * olarak �iz
                 
                    %Clusterhead olarak se�ilen node ile baz istasyonunun birbirine olan uzakl�g� hesapland� ve distance degi�kenine atand�
                    distance = sqrt((S(i).xd - (S(n + 1).xd)) ^ 2 + (S(i).yd - (S(n + 1).yd)) ^ 2); %�klit ile mesafe hesab�
                    C(cluster).distance = distance; %Cluster �n distance kolonuna uzakl�g� atad�k.
                    C(cluster).id = i; %Cluster�n id sine node numaras�n� atad�k
                    X(cluster) = S(i).xd; %Cluster�n x kordinat�n� X dizisine ata
                    Y(cluster) = S(i).yd; %Cluster�n y kordinat�n� Y dizisine ata
                    cluster = cluster + 1; %Cluster� 1 art�rd�
                 
                    %Sens�r d���mleri k-bit verilerini iletti�inde, enerji da��t�m� ��yledir:
                 
                    % ClusterHeadin baz istasyonuna olan mesafesi ve do de�eri
                    % dikkate al�narak ne kadar enerji harcamas� gerektigine karar
                    % veriliyor ve  ch in Enerjisni if sat�r�na g�re ve form�ldeki �ekilde g�ncelledik
                    distance;
                    if (distance > do)
                        S(i).E = S(i).E - ((ETX + EDA) * (k) + Emp * k * (distance * distance * distance * distance));
                        S(i).E = S(i).E - ((ETX + EDA) * (k) + Emp * k * (distance * distance * distance * distance));
                    end
                    if (distance <= do)
                        S(i).E = S(i).E - ((ETX + EDA) * (k) + Efs * k * (distance * distance));
                        S(i).E = S(i).E - ((ETX + EDA) * (k) + Efs * k * (distance * distance));
                    end
                end
             
            end
        end
    end
 
    % Burada STATISTICS ve CLUSTERHS dizisine �uanki rounda olu�turulan clusterHead say�s�n� yazd�k
    % (1 eksiltmemizin sebebi cluster de�i�keninin 1 den ba�lamas�, dizi indeksi 1 den ba�lad��� i�in)
    STATISTICS(r + 1).CLUSTERHEADS = cluster - 1;
    CLUSTERHS(r + 1) = cluster - 1;
 
    %Normal nodlar�n hangi cluster head e baglanacag�n� se�iyoruz.
 
    for i = 1:1:n
        if (S(i).type == 'N' && S(i).E > 0)
         
            %A�ag�daki sat�rlarda for d�ng�s�nde baz istanyonu da dahil olmak �zere t�m clusterheadlerin elimizdeki noda olan uzakl�g�na bakarak mesafe ve node numaras�n� tespit ettik
            %ve min_dis degi�keninde mesafayi, min_dis_cluster de�i�keninde ise bu clusterheadin numaras�n� tuttuk
            if (cluster - 1 >= 1) %Cluster say�s� birden fazlam�. cluster-1 ise ger�ek cluster say�s� normalde 0 tan�mlanamad�g� i�in
                min_dis = sqrt((S(i).xd - S(n + 1).xd) ^ 2 + (S(i).yd - S(n + 1).yd) ^ 2); %Elimizdeki node un baz istasyonuna olan uzakl�g� hesapalanarak mindis degi�kenine atand�
                min_dis_cluster = 1; %Degi�kenken olu�turmu� 1 den ba�latm��
                for c = 1:1:cluster - 1 %Clusterhead kadar d�n
                    temp = min(min_dis, sqrt((S(i).xd - C(c).xd) ^ 2 + (S(i).yd - C(c).yd) ^ 2)); %Burda node kendisine baz istasyonu mu yak�n yoksa var olan clusterheadlarden biri mi yak�n onu sorguluyo
                    if (temp < min_dis)
                        min_dis = temp;
                        min_dis_cluster = c; %Bize en yak�n olan clusterhead in numaras� minimum cluster degi�keninde tut.
                    end
                end
             
                %Elimizdeki node az �nce tespit etti�imiz kendisine en yak�n Clusterheade paket g�ndererek enerji harc�yor
                min_dis; %min_dis degi�keni do degerinden b�y�kse 1.if deki kadar enerji harc�yor
                if (min_dis > do)
                    S(i).E = S(i).E - (ETX * (k) + Emp * k * (min_dis * min_dis * min_dis * min_dis));
                end
                if (min_dis <= do); %min_dis degi�keni do degerinden k�y�kse 2.if deki kadar enerji harc�yor
                    S(i).E = S(i).E - (ETX * (k) + Efs * k * (min_dis * min_dis));
                end
             
                %Cluster headlerin baz veri g�ndererek enerjilerinin g�ncellenmesi
                if (min_dis > 0)
                    %Elimizdeki clusterhead ile baz istasyonu aras�ndaki mesafeyi hesaplad�k
                    distance = sqrt((S(C(min_dis_cluster).id).xd - (S(n + 1).xd)) ^ 2 + (S(C(min_dis_cluster).id).yd - (S(n + 1).yd)) ^ 2);
                 
                    S(C(min_dis_cluster).id).E = S(C(min_dis_cluster).id).E - ((ERX + EDA) * k); %ClusterHeadin yukarda ki node dan veri toplarken kaybetti�i enerjiyi g�ncelledik
                 
                    %Clusterhead baz istasyonuna olan uzakl�g�na g�re baz istanyonuna veri g�nderip enerji harc�yor
                    if (distance > do)
                        S(C(min_dis_cluster).id).E = S(C(min_dis_cluster).id).E - ((ETX + EDA) * (k) + Emp * k * (distance * distance * distance * distance));
                    end
                    if (distance <= do)
                        S(C(min_dis_cluster).id).E = S(C(min_dis_cluster).id).E - ((ETX + EDA) * (k) + Efs * k * (distance * distance));
                    end
                 
                    PACKETS_TO_CH(r + 1) = n - total_dead - cluster + 1; %?
                end
             
                S(i).min_dis = min_dis; % node dizisinin,i�eresinde elimizdeki nodun min_dis �zelli�i�e kendisine en yak�n olan cluster�nhead in mesafesini atat�k
                S(i).min_dis_cluster = min_dis_cluster; % node dizisinin,i�eresinde elimizdeki nodun min_dis_cluster �zelli�i�e kendisine en yak�n olan cluster�nhead in id degerini atat�k
             
            end
        end
    end
    hold on;
 
    countCHs;
    rcountCHs = rcountCHs + countCHs;
 
    %Her bir roun i�in, elimizdeki nodelar�n enerjileri toplay�p node say�s�na b�lerek o roundun ortalama enerciyi hesaplad�k
    sum = 0;
    for i = 1:1:n
        if (S(i).E > 0)
            sum = sum + S(i).E;
        end
    end
    avg = sum / n;
    STATISTICS(r + 1).AVG = avg; %STATISTICS degi�keninde rounda ait ortalama enerjiyi yazd�k
    sum;
 
    %
    % Code for Voronoi Cells
    % Unfortynately if there is a small
    % number of cells, Matlab's voronoi
    % procedure has some problems
    %
    % [vx,vy]=voronoi(X,Y);
    % plot(X,Y,'r*',vx,vy,'b-');
    % hold on;
    %
    % voronoi(X,Y);
    % axis([0 xm 0 ym]);
 
end
plot(S(n + 1).xd, S(n + 1).yd, 'x'); %Baz istasyonunu ekrana x olarak �izdirdik

figure(3);
for r = 0:1:rmax - 1
    ylabel('Her Nodeun Ortalama Enerjisi');
    xlabel('Round Numaras�');
    plot([r r + 1], [STATISTICS(r + 1).AVG STATISTICS(r + 2).AVG], 'red');
    hold on;
end

figure(4);
for r = 0:1:rmax - 1
    ylabel('�l� Node Say�s�');
    xlabel('Round Numaras�');
    plot([r r + 1], [STATISTICS(r + 1).DEAD STATISTICS(r + 2).DEAD], 'red');
    hold on;
end